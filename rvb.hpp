#ifndef RVB_HPP
#define RVB_HPP

//rvb.cpp
//Class that build the rvb state as a dimer + spin

#include <vector>
#include <iostream>
#include <random>

//#include "hypercube.hpp"
#include "square_lattice.hpp"

template<class Lattice> class RVB{

private:
  
  Lattice &lattice_;
  std::mt19937 rgen_;

public:
  int numSpins_; //total number of sites
  int D_;
  //int numVB_;
  int L_;
  int num_loops_;

  std::vector<int> spins_;
  std::vector<int> VB_;
  
  RVB(Lattice &lattice,int seed):lattice_(lattice){
    Init(seed); 
  }
  
  inline Lattice GetLattice(){return lattice_;}

  void Init(int seed){
    L_ = lattice_.LinSize();
    numSpins_  = lattice_.Nsites();
    //numVB_ = numSpins_/2;
    spins_.resize(numSpins_);
    VB_.resize(numSpins_);
    rgen_.seed(seed);
    //rgen_.seed();
    
    // Reset the configuration into the trivial topological sector
    for (int i=0;i<numSpins_;i+=2){
      VB_[i]  = lattice_.neighbours_[i][0];
      VB_[lattice_.neighbours_[i][0]]= i;
    }
    std::uniform_real_distribution<double> distribution(0,1);
    int spin_val;
    for (int i=0;i<numSpins_;i+=2){
      spin_val = (distribution(rgen_) < 0.5);
      spins_[i] = spin_val;
      spins_[VB_[i]] = spin_val ^ 1;
    }
  }

  void CopyVBToAncilla(){
    for(int i=0;i<VB_.size()/2;i++){
      VB_[i+VB_.size()/2] = VB_[i]+VB_.size()/2;
    }
  }
  void SetSpins(std::vector<int> spins){
    for(int i=0;i<spins.size();i++){
      spins_[i] = spins[i];
    }
  }
  void SetVB(std::vector<int> dimers){
    for(int i=0;i<VB_.size();i++){
      VB_[i] = VB_[i];
    }
  }
  void SetState(RVB &rvb){
    SetSpins(rvb.spins_);
    SetVB(rvb.VB_);
  }

  int LocalBondUpdate(int plaq){
    int flag = 0;
    //Check if two horizontal VBs in the plaquette
    if(VB_[lattice_.SitesOnPlaquettes_[plaq][0]] == lattice_.SitesOnPlaquettes_[plaq][1] &&
       VB_[lattice_.SitesOnPlaquettes_[plaq][2]] == lattice_.SitesOnPlaquettes_[plaq][3]){
      //Check if spins configurations are compatible
      if(spins_[lattice_.SitesOnPlaquettes_[plaq][0]] != spins_[lattice_.SitesOnPlaquettes_[plaq][3]] &&
         spins_[lattice_.SitesOnPlaquettes_[plaq][1]] != spins_[lattice_.SitesOnPlaquettes_[plaq][2]]) {
        // Perform the update
        //std::cout<<"Horizontal Bond Loop Update ACCEPTED"<<std::endl;
        flag = 1;
        VB_[lattice_.SitesOnPlaquettes_[plaq][0]] = lattice_.SitesOnPlaquettes_[plaq][3];
        VB_[lattice_.SitesOnPlaquettes_[plaq][1]] = lattice_.SitesOnPlaquettes_[plaq][2];
        VB_[lattice_.SitesOnPlaquettes_[plaq][2]] = lattice_.SitesOnPlaquettes_[plaq][1];
        VB_[lattice_.SitesOnPlaquettes_[plaq][3]] = lattice_.SitesOnPlaquettes_[plaq][0];
      }
    }
    //Check if two vertical VBs in the plaquette
    else if(VB_[lattice_.SitesOnPlaquettes_[plaq][0]] == lattice_.SitesOnPlaquettes_[plaq][3] &&
            VB_[lattice_.SitesOnPlaquettes_[plaq][1]] == lattice_.SitesOnPlaquettes_[plaq][2]){
      //Check if spins configurations are compatible
      if(spins_[lattice_.SitesOnPlaquettes_[plaq][0]] != spins_[lattice_.SitesOnPlaquettes_[plaq][1]] &&
         spins_[lattice_.SitesOnPlaquettes_[plaq][2]] != spins_[lattice_.SitesOnPlaquettes_[plaq][3]]) {
        // Perform the update
        flag = 1;
        //std::cout<<"Vertical Bond Loop Update ACCEPTED"<<std::endl;
        VB_[lattice_.SitesOnPlaquettes_[plaq][0]] = lattice_.SitesOnPlaquettes_[plaq][1];
        VB_[lattice_.SitesOnPlaquettes_[plaq][1]] = lattice_.SitesOnPlaquettes_[plaq][0];
        VB_[lattice_.SitesOnPlaquettes_[plaq][2]] = lattice_.SitesOnPlaquettes_[plaq][3];
        VB_[lattice_.SitesOnPlaquettes_[plaq][3]] = lattice_.SitesOnPlaquettes_[plaq][2];
      }
    }
    return flag;
  }

  //Spin update compatible with beta
  void SpinUpdate(RVB &rvb,int &ratio){
    std::vector<int> already_updated;
    already_updated.assign(numSpins_,0);
    int next_site,spin_val;
    int num_loops = 0;
    std::uniform_int_distribution<int> distribution(0,1); 
    
    if (ratio) Swap();  //swap the states

    for(int i=0;i<numSpins_;i++){
      if (already_updated[i] == 0) {
        //std::cout<<"Starting spin = S(" << i <<") = "<< spins_[i]<< std::endl;
        spin_val = distribution(rgen_);
        spins_[i] = spin_val;
        //std::cout<<"  New value = S(" << i <<") = "<< spins_[i]<< std::endl;
        already_updated[i] = 1;
        next_site = VB_[i];
        //std::cout<<"Next spin = S(" << next_site<<") = "<< spins_[next_site] << std::endl;
        while(already_updated[next_site] == 0){
          if(already_updated[next_site] == 1) 
            std::cout<<"Loop Error 1"<<std::endl;
          else 
            already_updated[next_site] = 1;
          spin_val = spin_val ^ 1;
          spins_[next_site] = spin_val;
          //std::cout<<"  New value = S(" << next_site<<") = "<< spins_[next_site] << std::endl;

          next_site = rvb.VB_[next_site];
          //std::cout<<"Next spin = S(" << next_site<<") = "<< spins_[next_site] << std::endl;
          if(already_updated[next_site] != 1)
            already_updated[next_site] = 1;
          else break;

          spin_val = spin_val ^ 1;
          spins_[next_site] = spin_val;
          next_site = VB_[next_site];
          //std::cout<<"Next spin = " << next_site << std::endl;
        }
        num_loops++;
      }
    }
    //Set the spin state of rvb
    rvb.SetSpins(spins_);
    //Save the number of loops
    num_loops_ = num_loops;
    rvb.num_loops_ = num_loops;
    
    if (ratio) Swap();  //unswap the states
  }

  int Overlap(RVB &rvb) {
    std::vector<int> already_counted;
    already_counted.assign(numSpins_,0);
    int next_site;
    int num_loops = 0;

    for(int i=0;i<numSpins_;i++){
      if (already_counted[i] == 0) {
        already_counted[i] = 1;
        next_site = VB_[i];
        while(already_counted[next_site] == 0){
          if(already_counted[next_site] == 1) 
            std::cout<<"Loop Error 1"<<std::endl;
          else 
            already_counted[next_site] = 1;
          next_site = rvb.VB_[next_site];
          //std::cout<<"Next spin = S(" << next_site<<") = "<< spins_[next_site] << std::endl;
          if(already_counted[next_site] != 1)
            already_counted[next_site] = 1;
          else break;
          next_site = VB_[next_site];
        }
        num_loops++;
      }
    }
    return num_loops;
  }

  void Swap(std::vector<int> regionA) {
    
    int site,replica;
    int inBondWithSite,inBondWithReplica;
    int tmp;

    for(int i=0;i<numSpins_/2;i++){
      
      if(regionA[i] == 1){
        site = i;
        replica = numSpins_/2 + i;
        // Swap the spin state
        tmp = spins_[site];
        spins_[site] = spins_[replica];
        spins_[replica] = tmp;

        //Swap the bond state
        inBondWithSite    = VB_[site];
        inBondWithReplica = VB_[replica];
      
        VB_[site]    = inBondWithReplica;
        VB_[replica] = inBondWithSite;
        VB_[inBondWithSite] = replica;
        VB_[inBondWithReplica] = site;

        //if(regionA[bondB] == 1){
        //  VB_[siteA] = bondB - numSpins_/2;
        //}
        //else{
        //  VB_[siteA] = bondB;
        //  VB_[bondB] = siteA;
        //}
        //if(regionA[bondA] == 1){
        //  VB_[siteB] = bondA + numSpins_/2;
        //}
        //else{
        //  VB_[siteB] = bondA;
        //  VB_[bondA] = siteB;
        //}
      }
    }
  }

  int TopologicalSectorRealX(){
    int sector = 0;
    for(int i=0;i<L_;i+=2){
      if(VB_[i] == lattice_.neighbours_[i][1])
        sector += 1;
      if (VB_[i+1] == lattice_.neighbours_[i+1][1])
        sector -= 1;
    }
    return sector;
  }

  int TopologicalSectorRealY(){
    int sector = 0;
    for(int i=0;i<L_;i+=2){
      if(VB_[L_*i] == lattice_.neighbours_[L_*i][0])
        sector += 1;
      if (VB_[L_*(i+1)] == lattice_.neighbours_[L_*(i+1)][0])
        sector -= 1;
    }
    return sector;
  }

  int TopologicalSectorReplicaX(){
    int sector = 0;
    for(int i=0;i<L_;i+=2){
      if(VB_[numSpins_/2+i] == lattice_.neighbours_[numSpins_/2+i][1])
        sector += 1;
      if (VB_[numSpins_/2+i+1] == lattice_.neighbours_[numSpins_/2+i+1][1])
        sector -= 1;
    }
    return sector;
  }

  int TopologicalSectorReplicaY(){
    int sector = 0;
    for(int i=0;i<L_;i+=2){
      if(VB_[numSpins_/2+L_*i] == lattice_.neighbours_[numSpins_/2+L_*i][0])
        sector += 1;
      if (VB_[numSpins_/2+L_*(i+1)] == lattice_.neighbours_[numSpins_/2+L_*(i+1)][0])
        sector -= 1;
    }
    return sector;
  }
 
  int CheckTopologicalSectors(int &Wx, int &Wy){
    if (TopologicalSectorRealX() != Wx) return 1;
    else if (TopologicalSectorRealY() != Wy) return 1;
    else if (TopologicalSectorReplicaX() != Wx) return 1;
    else if (TopologicalSectorReplicaY() != Wy) return 1;
    else return 0;
  }

  void print(){

    std::cout<<"VB basis: \n";
    for (int i=0;  i<VB_.size(); i++)
        //VBasis[i].print();
        std::cout<<i<<"->"<<VB_[i]<<std::endl;
  };//print


  void printTOPO(){ //print the topological sectors
    std::cout<<"("<<TopologicalSectorRealX()<<","<<TopologicalSectorRealY()<<")"<<std::endl;
  }//printTOPO



    //void PrintSpins(){
    //  for(int y=0;y< lattice_.LinSize(); y++){
    //    for(int x=0;x< lattice_.LinSize(); x++){
    //      if(spins_[lattice_.Index(x,y)] == 1)
    //        PRINT_RED("+");
    //      else
    //        PRINT_GREEN("-");
    //      std::cout<<"  ";
    //    }
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl<<std::endl;
    //}

    //void Print(){
    //  for(int y=0;y< lattice_.LinSize(); y++){
    //    for(int x=0;x< lattice_.LinSize(); x++){
    //      if(spins_[lattice_.Index(x,y)] == 1)
    //        PRINT_RED("o");
    //      else
    //        PRINT_GREEN("o");

    //      if(dimers_[lattice_.LinksOnSites_[lattice_.Index(x,y)][0]]==1){
    //        PRINT_BLUE("---");
    //      }
    //      else{
    //        std::cout<<"   ";
    //      }
    //    }
    //    std::cout<<std::endl;
    //    for(int x=0;x< lattice_.LinSize(); x++){
    //      if(dimers_[lattice_.LinksOnSites_[lattice_.Index(x,y)][1]]==1){
    //        PRINT_BLUE("|   ");
    //      }
    //      else{
    //        std::cout<<"    ";
    //      }
    //    }
    //    std::cout<<std::endl;
 
    //  }
    //  std::cout<<std::endl<<std::endl;
    //}

    //
    //  
    //void TestBondUpdate(int plaq){
    //
    //  std::cout<<" BOND UPDATE: attempt on plaquette "<< plaq;
    //  std::vector<int> links;
    //  for(int i=0; i<4; i++){
    //    if(dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] == 1)
    //      links.push_back(lattice_.LinksOnPlaquettes_[plaq][i]);
    //  }
    //  
    //  //std::cout<<"Number of dimers = " << tmp << "   ";
    //  if ((links.size() == 2) && (spins_[lattice_.SitesOnLinks_[links[0]][0]] != spins_[lattice_.SitesOnLinks_[links[1]][0]])){
    //    for (int i=0; i<4; i++){
    //      dimers_[lattice_.LinksOnPlaquettes_[plaq][i]] = 1 - dimers_[lattice_.LinksOnPlaquettes_[plaq][i]];
    //    }
    //    std::cout<< " ACCEPTED" << std::endl;
    //    Print();
    //  }
    //  else std::cout<< " REJECTED" << std::endl;
    //}
 
};

#endif
