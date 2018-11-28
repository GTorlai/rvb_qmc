#ifndef QMC_HPP
#define QMC_HPP

#include <vector>
#include <iostream>
#include <random>
#include <math.h> 
#include "rvb.hpp"
#include "parameters.hpp"

template<class Lattice> class QMC{

private:
  
  Lattice &lattice_;
  RVB<Lattice> bra_;
  RVB<Lattice> ket_;
  std::mt19937 rgen_;

public:
  std::vector<std::vector<double> > SpinSpinCorrelation_;


  int numSpins_; //total number of sites
  int Ndimers_;
  int nsamples_node_;
  int nburn_;
  //int flag_;
  int Wx_,Wy_;
 
  




  //Functions
  QMC(Lattice &lattice,Parameters &pars):lattice_(lattice),
                        bra_(lattice,pars.seed_bra_),
                        ket_(lattice,pars.seed_ket_),
                        nburn_(pars.nburn_){
    numSpins_ = lattice_.Nsites();
    Wx_ = pars.Wx_;
    Wy_ = pars.Wy_;
    rgen_.seed(pars.seed_qmc_);
    nsamples_node_ = std::ceil(double(pars.nMC_));// / double(pars.totalnodes_));
    Init();
  }

  void Init(){
    ket_.SpinUpdate(bra_);
    ket_.CopyVBToAncilla();
    bra_.SetState(ket_);
    std::cout<<"Ket topological sector: ";
    ket_.printTOPO();
    std::cout<<"Bra topological sector: ";
    bra_.printTOPO();
  }
 
  void QMCrun() {
    for (int i=0;i<nburn_;i++){
      Sweep();
    }
    for (int i=0;i<nsamples_node_;i++){
      //std::cout<<"Sweep # "<<i<<std::endl<<std::endl;
      Sweep();
      GetSpinSpinCorrelation();
    }
  }

  void Sweep(){
    //TODO 3d CHANGE
    std::uniform_int_distribution<int> dist(0,numSpins_-1);
    int plaq;

    //TODO 3d CHANGE
    for(int p=0;p<numSpins_;p++){
      plaq = dist(rgen_);
      ket_.LocalBondUpdate(plaq);
      plaq = dist(rgen_);
      bra_.LocalBondUpdate(plaq);
    }
    //ket_.print();
    ket_.SpinUpdate(bra_);
    //PrintTransitionGraph();
    SanityCheck();
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
