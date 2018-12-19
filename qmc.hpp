#ifndef QMC_HPP
#define QMC_HPP

#include <vector>
#include <iostream>
#include <random>
#include <math.h> 
#include <string>
#include "rvb.hpp"
#include "parameters.hpp"
#include "utils.hpp"
template<class Lattice> class QMC{

private:
  
  Lattice &lattice_;
  RVB<Lattice> bra_;
  RVB<Lattice> ket_;
  std::mt19937 rgen_;
  
  std::vector<int> regionA_;
  //std::vector<std::vector<int> > regionA_;
  std::vector<int> regionX_;
  int ratio_;

public:
  std::vector<std::vector<double> > SpinSpinCorrelation_;
  //std::vector<std::vector<double> > RenyiEntropy_;
  std::vector<double> RenyiEntropy_;

  int numSpins_; //total number of sites
  int numPlaqs_; //total number of plaquettes
  int Ndimers_;
  unsigned long long int nsamples_node_;
  int nburn_;
  int Wx_,Wy_;
  int nloops_den,nloops_num;

  //Functions
  QMC(Lattice &lattice,Parameters &pars):lattice_(lattice),
                        bra_(lattice,pars.seed_bra_),
                        ket_(lattice,pars.seed_ket_),
                        nburn_(pars.nburn_){
    numSpins_ = lattice_.Nsites();
    numPlaqs_ = lattice_.Nplaqs();
    Wx_ = pars.Wx_;
    Wy_ = pars.Wy_;
    ratio_ = pars.ratio_;
    rgen_.seed(pars.seed_qmc_);
    //nsamples_node_ = std::ceil(double(pars.nMC_)/double(pars.totalnodes_));
    nsamples_node_ = pars.nMC_;
    Init();
  }

  void Init(){
    ket_.SpinUpdate(bra_);
    ket_.CopyVBToAncilla();
    bra_.SetState(ket_);
  }
 
  void Reset(){
    regionA_.clear();
    regionX_.clear();
    SpinSpinCorrelation_.clear();
    RenyiEntropy_.clear();
  }

  void LoadRegion(Parameters &pars,int regionID){
    std::string fname;
    std::string path = SimulationName(pars) + "/regions/";
    fname = path + "regionA_" + std::to_string(regionID) + ".txt";
    std::ifstream fin(fname);
    int tmp;
    for(int i=0;i<numSpins_/2;i++){
      fin >> tmp;
      regionA_.push_back(tmp);
    }
    fin.close();
    
    if(ratio_){
      fname = path+"regionX_" + std::to_string(regionID) + ".txt";
      std::ifstream fin(fname);
      int tmp;
      for(int i=0;i<numSpins_/2;i++){
        fin >> tmp;
        regionX_.push_back(tmp);
      }
      fin.close();
      ket_.SetRegion(regionX_);
    }
  }

  void Equilibrate() {
    for (int i=0;i<nburn_;i++){
      Sweep();
    }
  }
  
  void Run() {
    for (int i=0;i<nsamples_node_;i++){
      Sweep();
      //GetSpinSpinCorrelation();
      GetEntanglementEntropy(); 
    }
  }

  void Sweep(){
    std::uniform_int_distribution<int> dist(0,numPlaqs_-1);
    int plaq;

    for(int p=0;p<numPlaqs_/2;p++){
      plaq = dist(rgen_);
      ket_.LocalBondUpdate(plaq);
      plaq = dist(rgen_);
      bra_.LocalBondUpdate(plaq);
    }
    //ket_.print();
    ket_.SpinUpdate(bra_);
    //PrintTransitionGraph();
    //SanityCheck();
    if (ket_.CheckTopologicalSectors(Wx_,Wy_) || bra_.CheckTopologicalSectors(Wx_,Wy_)){
      std::cout<<"WRONG TOPOLOGICAL SECTOR"<<std::endl;
      exit(0);
    }
  }

  bool CheckLoopSharing(int siteA,int siteB){
    int next_site,counter;
    bool flag = false;
    next_site = siteA;
    do {
      next_site = ket_.VB_[next_site];
      if(next_site == siteB){
        flag=true;
        break;
      }
      next_site = bra_.VB_[next_site];
      if(next_site == siteB){
        flag=true;
        break;
      }
    } while(next_site != siteA);
    return flag;
  }

  void GetSpinSpinCorrelation() {
    bool flag;
    std::vector<double> tmp(lattice_.LinSize()/2);
    for(int x=1;x<lattice_.LinSize()/2+1;x++){
      flag = CheckLoopSharing(0,x); 
      if (flag==false){
        tmp[x-1] = 0.0;
      }
      else{
        if (x % 2 == 0)
          tmp[x-1] = 0.75;
        else
          tmp[x-1] = -0.75;
      }
      //std::cout<<tmp[x-1] << "  ";
    }
    //std::cout<<std::endl;
    SpinSpinCorrelation_.push_back(tmp);
  }

  void GetEntanglementEntropy(){
    //Measure the swap operator
    if(ratio_){
      ket_.Swap(regionX_);
      nloops_den = ket_.Overlap(bra_);
      ket_.Swap(regionX_);
    }
    else{
      nloops_den = ket_.Overlap(bra_);
    }
    ket_.Swap(regionA_);
    nloops_num = ket_.Overlap(bra_);
    ket_.Swap(regionA_);
    RenyiEntropy_.push_back(1.0*std::pow(2,nloops_num-nloops_den));
  }









  ////--------- TEST ----------//


  void SanityCheck(){
    
    //Check that the two spin configuration matches
    bool flag = true;

    for(int i=0;i<numSpins_;i++){
      if (ket_.spins_[i] != bra_.spins_[i]){
        flag=false;
        std::cout<< "ERROR: bra and ket do not match" <<std::endl;
        exit(0);
      }
    }

    //Check that the dimers are compatible with spins
    flag=true;
    for(int i=0;i<ket_.VB_.size();i++){
      if(ket_.spins_[i] == ket_.spins_[ket_.VB_[i]]){
        flag=false;
        std::cout<< "ERROR: dimers and spins not compatible in the ket  " << i << std::endl;
        exit(0);
      }
      if(bra_.spins_[i] == bra_.spins_[bra_.VB_[i]]){
        flag=false;
        std::cout<< "ERROR: dimers and spins not compatible in the bra  " << i << std::endl;
        exit(0);
      }
    }
  }

  void PrintState(RVB<Lattice> &rvb){
    for(int y=0;y< lattice_.LinSize(); y++){
      for(int x=0;x< lattice_.LinSize(); x++){
        if(rvb.spins_[lattice_.Index(x,y)] == 1)
          PRINT_RED("o");
        else
          PRINT_GREEN("o");

        if(rvb.VB_[lattice_.Index(x,y)] == lattice_.Index(x+1,y)){
          PRINT_BLUE("---");
        }
        else{
          std::cout<<"   ";
        }
      }
      std::cout<<std::endl;
      for(int x=0;x< lattice_.LinSize(); x++){
        if(rvb.VB_[lattice_.Index(x,y)] == lattice_.Index(x,y+1)){
          PRINT_BLUE("|   ");
        }
        else{
          std::cout<<"    ";
        }
      }
      std::cout<<std::endl;
 
    }
    std::cout<<std::endl<<std::endl;
  }

  void PrintTransitionGraph(){
    for(int y=0;y< lattice_.LinSize(); y++){
      for(int x=0;x< lattice_.LinSize(); x++){
        if(ket_.spins_[lattice_.Index(x,y)] == 1)
          PRINT_RED("o");
        else
          PRINT_GREEN("o");

        if(ket_.VB_[lattice_.Index(x,y)] == lattice_.Index(x+1,y) && bra_.VB_[lattice_.Index(x,y)] == lattice_.Index(x+1,y)){
          PRINT_BLUE("===");
        }
        else if(ket_.VB_[lattice_.Index(x,y)] == lattice_.Index(x+1,y) || bra_.VB_[lattice_.Index(x,y)] == lattice_.Index(x+1,y)){
          PRINT_BLUE("---");
        }
        else{
          std::cout<<"   ";
        }
      }
      std::cout<<std::endl;
      for(int x=0;x< lattice_.LinSize(); x++){
        if(ket_.VB_[lattice_.Index(x,y)] == lattice_.Index(x,y+1) && bra_.VB_[lattice_.Index(x,y)] == lattice_.Index(x,y+1)){
          PRINT_BLUE("||  ");
        }
        else if(ket_.VB_[lattice_.Index(x,y)] == lattice_.Index(x,y+1) || bra_.VB_[lattice_.Index(x,y)] == lattice_.Index(x,y+1)){
          PRINT_BLUE("|   ");
        }
        
        else{
          std::cout<<"    ";
        }
      }
      std::cout<<std::endl;
 
    }
    std::cout<<std::endl<<std::endl;
  }


};

#endif
