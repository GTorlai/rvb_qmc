#ifndef QMC_HPP
#define QMC_HPP

//rvb.cpp
//Class that build the rvb state as a dimer + spin

#include <vector>
#include <iostream>
#include <random>
#include <math.h> 
#include "rvb.hpp"

template<class Lattice> class QMC{

private:
  
  Lattice &lattice_;
  RVB<Lattice> bra_;
  RVB<Lattice> ket_;
  std::mt19937 rgen_;

  std::vector<int> TransitionGraph_;
  
public:
  std::vector<double> SpinSpinCorrelation_;
  std::vector<double> DimerDimerCorrelation_;
  std::vector<double> SpinSpinCorrelation_Var_;
  std::vector<double> DimerDimerCorrelation_Var_;
  std::vector<double> SpinSpinCorrelation_Avg_;
  std::vector<double> DimerDimerCorrelation_Avg_;

  //int D_; //dimension
  int Nspins_; //total number of sites
  int Ndimers_;
  //Functions
  QMC(Lattice &lattice,int seedQMC,int seedBra,int seedKet):lattice_(lattice),
                        bra_(lattice,seedBra),
                        ket_(lattice,seedKet){
    Init(seedQMC); 
  }

  void Init(int seedQMC){
    Nspins_ = lattice_.Nsites();
    SpinSpinCorrelation_.resize(Nspins_-1);
    SpinSpinCorrelation_Avg_.resize(Nspins_-1);
    SpinSpinCorrelation_Var_.resize(Nspins_-1);
    std::fill(SpinSpinCorrelation_.begin(), SpinSpinCorrelation_.end(), 0.0);
    std::fill(SpinSpinCorrelation_Avg_.begin(), SpinSpinCorrelation_Avg_.end(), 0.0);
    std::fill(SpinSpinCorrelation_Var_.begin(), SpinSpinCorrelation_Var_.end(), 0.0);
    
    //Ndimers_ = int(lattice_.Nlinks()/2);
    //ket_.PrintSpins();
    //bra_.PrintSpins();
    
    TransitionGraph_.resize(lattice_.Nlinks());
    rgen_.seed(seedQMC);
    Reset();
  }
  
  void Reset(){
    ket_.Reset();
    bra_.SetSpins(ket_.spins_);
    bra_.SetDimers(ket_.dimers_);
    //bra_.Reset();
  }
   
  void GetTransitionGraph(){
    for(int i=0;i<lattice_.Nlinks();i++){ 
      //if (ket_.dimers_[i] == bra_.dimers_[i])
      //  TransitionGraph_[i] = 0;
      //else
      //  TransitionGraph_[i] = 1;
      TransitionGraph_[i] = ket_.dimers_[i] + bra_.dimers_[i];
    }
  }

  void SpinUpdate(){
    std::uniform_int_distribution<int> distributionINT(0,Nspins_);
    std::uniform_real_distribution<double> distributionREAL(0,1);
    int site,site0,counter;
    int next_site;
    
    for(int l=0;l<lattice_.Nlinks();l++){
      if (TransitionGraph_[l] != 0) {
        site0 = lattice_.SitesOnLinks_[l][0]; 
        site=site0;
        next_site = -1;
        if (distributionREAL(rgen_) > 0.5){
          do {
            counter = 0;
            ket_.spins_[site] = -ket_.spins_[site];
            for(int i=0;i<lattice_.Neighbours_.size();i++){
              next_site = lattice_.Neighbours_[site][i];
              if(TransitionGraph_[lattice_.LinksOnSites_[site][i]] != 0){
                counter = i;
                break;
              }
            }
            TransitionGraph_[lattice_.LinksOnSites_[site][counter]]--;
            site = next_site;
          } while(site != site0);
        }
      }
    }
  }
 
  void QMCrun(int nburn,int nsweeps) {
    for (int i=0;i<nburn;i++){
      Sweep();
    }
    for (int i=0;i<nsweeps;i++){
      Sweep();
      GetTransitionGraph();
      GetSpinSpinCorrelation();
      UpdateMeasurement(i+1,SpinSpinCorrelation_Avg_,SpinSpinCorrelation_Var_,SpinSpinCorrelation_);
    }
  }

  void Sweep(){
    //TODO 3d CHANGE
    std::uniform_int_distribution<int> dist(0,Nspins_-1);
    int plaq;

    //TODO 3d CHANGE
    //bra_.Print();
    for(int p=0;p<Nspins_;p++){
      plaq = dist(rgen_);
      ket_.BondUpdate(plaq);
      plaq = dist(rgen_);
      bra_.BondUpdate(plaq);
    }
    GetTransitionGraph();
    //PrintTransitionGraph();
    SpinUpdate();
    bra_.SetSpins(ket_.spins_);
    SanityCheck();
      
  }

  
  void UpdateMeasurement(int count,std::vector<double> &avg,std::vector<double> &M2,std::vector<double> &newvalue){
    //count += 1;
    //for (int i=0;i<avg.size();i++){
    //  avg[i] += newvalue[i];
    //}
    std::vector<double> delta,delta2;
    delta.resize(avg.size());
    delta2.resize(avg.size());
    
    for (int i=0;i<avg.size();i++){
      delta[i] = newvalue[i] - avg[i];
      avg[i] += delta[i] / float(count);
      delta2[i] = newvalue[i] - avg[i];
      M2[i] += delta[i] * delta2[i];
    }
  }

  void GetSpinSpinCorrelation() {
    bool flag;
    for(int x=1;x<lattice_.LinSize();x++){
      flag = CheckLoopSharing(8,8+x); 
      if (flag==false){
        SpinSpinCorrelation_[x-1] = 0.0;
      }
      else{
        if (x % 2 == 0)
          SpinSpinCorrelation_[x-1] = 0.75;
        else
          SpinSpinCorrelation_[x-1] = -0.75;
      }
    }
  }

  bool CheckLoopSharing(int siteA,int siteB){
    std::vector<int> TransitionTMP;
    TransitionTMP.resize(lattice_.Nlinks());
    int site,next_site,counter;
    bool flag = false;
    site=siteA;
    next_site = -1;
    TransitionTMP = TransitionGraph_;
    //std::cout<<"Starting site: " << siteA << std::endl;
    //std::cout<<"Final site: " << siteB << std::endl;
    do {
      counter = 0;
      for(int i=0;i<lattice_.Neighbours_.size();i++){
        next_site = lattice_.Neighbours_[site][i];
        if(TransitionTMP[lattice_.LinksOnSites_[site][i]] != 0){
          counter = i;
          break;
        }
      }
      TransitionTMP[lattice_.LinksOnSites_[site][counter]]--;
      site = next_site;
      //std::cout<<"  next site = " << site << std::endl;
      if(site == siteB){
        //std::cout<< "SUCCESS" << std::endl << std::endl;
        flag = true;
        break;
      }
    } while(site != siteA);
    return flag;
  }





  //--------- TEST ----------//


  void SanityCheck(){
    
    //Check that the two spin configuration matches
    bool flag = true;

    for(int i=0;i<Nspins_;i++){
      if (ket_.spins_[i] != bra_.spins_[i]){
        flag=false;
        std::cout<< "ERROR: bra and ket do not match" <<std::endl;
        exit(0);
      }
    }

    //Check that the dimers are compatible with spins
    flag=true;
    for(int l=0;l<lattice_.Nlinks();l++){ 
      if(ket_.dimers_[l] == 1){
        if (ket_.spins_[lattice_.SitesOnLinks_[l][0]] == ket_.spins_[lattice_.SitesOnLinks_[l][1]]){
          flag=false;
          std::cout<< "ERROR: dimers and spins not compatible in the ket  " << l << std::endl;
          exit(0);
        }
      }
      if(bra_.dimers_[l] == 1){
        if (bra_.spins_[lattice_.SitesOnLinks_[l][0]] == bra_.spins_[lattice_.SitesOnLinks_[l][1]]){
          flag=false;
          std::cout<< "ERROR: dimers and spins not compatible in the bra  " << l <<std::endl;
          bra_.Print();
          exit(0);
        }
      }
    }
  }

  void TestSpinUpdate(){
    std::uniform_int_distribution<int> distributionINT(0,Nspins_);
    std::uniform_real_distribution<double> distributionREAL(0,1);
    int site,site0,counter;
    int next_site;
    //std::vector<int> transit;
    //transit = TransitionGraph_;
    
    for(int l=0;l<lattice_.Nlinks();l++){
      if (TransitionGraph_[l] != 0) {
        std::cout<<"First link = " << l;
        site0 = lattice_.SitesOnLinks_[l][0]; 
        site=site0;
        next_site = -1;
        std::cout<<"  -  Starting site = "<<site0<<std::endl;
        //TransitionGraph_[l] = 0;
        if (distributionREAL(rgen_) > 0.5){
        //if (1>0){
          do {
            counter = 0;
            ket_.spins_[site] = -ket_.spins_[site];
            for(int i=0;i<lattice_.Neighbours_.size();i++){
              next_site = lattice_.Neighbours_[site][i];
              std::cout<<"Proposed site: "<<next_site;
              std::cout<<"  -  Link in that direction = " << lattice_.LinksOnSites_[site][i];
              std::cout<<"  -  Graph = " << TransitionGraph_[lattice_.LinksOnSites_[site][i]];
              std::cout<<std::endl;
              if(TransitionGraph_[lattice_.LinksOnSites_[site][i]] != 0){
                std::cout<<"NEXT SITE = " << next_site << std::endl;
                counter = i;
                break;
              }
            }
            TransitionGraph_[lattice_.LinksOnSites_[site][counter]]--;
            site = next_site;
            PrintTransitionGraph();
          } while(site != site0);
        }
        //PrintTransitionGraph();
      }
    }
  }

  
  void Test(int nsweeps){
    //TODO 3d CHANGE
    std::uniform_int_distribution<int> distributionINT(0,Nspins_-1);
    int plaq;
    
    for (int i=0;i<nsweeps;i++){
      //TODO 3d CHANGE
      for(int p=0;p<Nspins_;p++){
        plaq = distributionINT(rgen_);
        ket_.TestBondUpdate(plaq);
        plaq = distributionINT(rgen_);
        bra_.TestBondUpdate(plaq);
      }
      GetTransitionGraph();
      ket_.Print();
      bra_.Print();
      PrintTransitionGraph();
      TestSpinUpdate();
      bra_.SetSpins(ket_.spins_);
      SanityCheck();
      //ket_.PrintDimers();
    }
  }


  void PrintTransitionGraph(){
    for(int y=0;y< lattice_.LinSize(); y++){
      for(int x=0;x< lattice_.LinSize(); x++){
        if(ket_.spins_[lattice_.Index(x,y)] == 1)
          PRINT_RED("o");
        else
          PRINT_GREEN("o");

        if(TransitionGraph_[lattice_.LinksOnSites_[lattice_.Index(x,y)][0]]==1){
          PRINT_BLUE("---");
        }
        else if(TransitionGraph_[lattice_.LinksOnSites_[lattice_.Index(x,y)][0]]==2){
          PRINT_BLUE("===");
        }
        else{
          std::cout<<"   ";
        }
      }
      std::cout<<std::endl;
      for(int x=0;x< lattice_.LinSize(); x++){
        if(TransitionGraph_[lattice_.LinksOnSites_[lattice_.Index(x,y)][1]]==1){
          PRINT_BLUE("|   ");
        }
        else if(TransitionGraph_[lattice_.LinksOnSites_[lattice_.Index(x,y)][1]]==2){
          PRINT_BLUE("||  ");
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
