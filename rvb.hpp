#ifndef RVB_HPP
#define RVB_HPP

//rvb.cpp
//Class that build the rvb state as a dimer + spin

#include <vector>
#include <iostream>
#include <random>

//#include "hypercube.hpp"
#include "square_lattice.hpp"
#define PRINT_RED(x) std::cout << "\033[1;31m" << x << "\033[0m"   //<< " "
#define PRINT_BLUE(x) std::cout << "\033[1;34m" << x << "\033[0m"  //<< " "
#define PRINT_GREEN(x) std::cout << "\033[1;32m" << x << "\033[0m" //<< " "


template<class Lattice> class RVB{

private:
  
  Lattice &lattice_;
  std::mt19937 rgen_;
public:
    int Nspins_; //total number of sites
    int Ndimers_;

    std::vector<int> spins_;
    std::vector<int> dimers_;
    
    RVB(Lattice &lattice,int seed):lattice_(lattice){
      Init(seed); 
    }
    
    inline Lattice GetLattice(){return lattice_;}

    void Init(int seed){
      Nspins_ = lattice_.Nsites();
      spins_.resize(Nspins_);
      dimers_.resize(lattice_.Nlinks());
      
      rgen_.seed(seed);
      //rgen_.seed();
      Reset();
    }
    
    void SetSpins(std::vector<int> spins){
      spins_ = spins;
    }
    void SetDimers(std::vector<int> dimers){
      dimers_ = dimers;
    }
    void Reset(){
      // Reset the configuration into the W=(0,0) topological sector
      for (int i=0;i<lattice_.Nlinks();i++){
        if(i % 4 == 0) 
          dimers_[i] = 1;
        else
          dimers_[i] = 0;
      }

      std::uniform_real_distribution<double> distribution(0,1);
      // Reset the configuration of the spins at random
      for (int i=0;i<lattice_.Nlinks();i++){
        if (dimers_[i] == 1){
          int svalue = 1-2*(distribution(rgen_) < 0.5);
          spins_[lattice_.SitesOnLinks_[i][0]] = svalue;
          spins_[lattice_.SitesOnLinks_[i][1]] = -svalue;
        }
      }
    }
    
    void BondUpdate(int plaq){
    
      std::vector<int> links;
      for(int i=0; i<4; i++){
        if(dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] == 1)
          links.push_back(lattice_.LinksOnPlaquettes_[plaq][i]);
      }
      if ((links.size() == 2) && (spins_[lattice_.SitesOnLinks_[links[0]][0]] != spins_[lattice_.SitesOnLinks_[links[1]][0]])){
        for (int i=0; i<4; i++){
          dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] = 1 - dimers_[lattice_.LinksOnPlaquettes_[plaq][i]];
        }
      }
    }
    
    void PrintSpins(){
      for(int y=0;y< lattice_.LinSize(); y++){
        for(int x=0;x< lattice_.LinSize(); x++){
          if(spins_[lattice_.Index(x,y)] == 1)
            PRINT_RED("+");
          else
            PRINT_GREEN("-");
          std::cout<<"  ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl<<std::endl;
    }

    void Print(){
      for(int y=0;y< lattice_.LinSize(); y++){
        for(int x=0;x< lattice_.LinSize(); x++){
          if(spins_[lattice_.Index(x,y)] == 1)
            PRINT_RED("o");
          else
            PRINT_GREEN("o");

          if(dimers_[lattice_.LinksOnSites_[lattice_.Index(x,y)][0]]==1){
            PRINT_BLUE("---");
          }
          else{
            std::cout<<"   ";
          }
        }
        std::cout<<std::endl;
        for(int x=0;x< lattice_.LinSize(); x++){
          if(dimers_[lattice_.LinksOnSites_[lattice_.Index(x,y)][1]]==1){
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

    
      
    void TestBondUpdate(int plaq){
    
      std::cout<<" BOND UPDATE: attempt on plaquette "<< plaq;
      std::vector<int> links;
      for(int i=0; i<4; i++){
        if(dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] == 1)
          links.push_back(lattice_.LinksOnPlaquettes_[plaq][i]);
      }
      
      //std::cout<<"Number of dimers = " << tmp << "   ";
      if ((links.size() == 2) && (spins_[lattice_.SitesOnLinks_[links[0]][0]] != spins_[lattice_.SitesOnLinks_[links[1]][0]])){
        for (int i=0; i<4; i++){
          dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] = 1 - dimers_[lattice_.LinksOnPlaquettes_[plaq][i]];
        }
        std::cout<< " ACCEPTED" << std::endl;
        Print();
      }
      else std::cout<< " REJECTED" << std::endl;
    }
 
};

#endif
