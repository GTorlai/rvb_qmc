#ifndef RVB_HPP
#define RVB_HPP

//rvb.cpp
//Class that build the rvb state as a dimer + spin

#include <vector>
#include <iostream>
#include <random>

#include "cubic_lattice.hpp"
#include "square_lattice.hpp"

template<class Lattice> class RVB{

private:
  
  Lattice &lattice_;          // Lattice object
  std::mt19937 rgen_;         // Random nunber generator

public:
  int numSpins_;              // Number of spins
  int D_;                     // Lattice dimension
  int L_;                     // Linear size of the lattice
  int num_loops_;             // Number of loops in the transition graph
  int ratio_;                 // Ratio flag
  std::vector<int> regionX_;  // RegionX for the ratio trick

  std::vector<int> spins_;    // Vector of spin values
  std::vector<int> VB_;       // Vector of valence bonds
  
  // Constructor
  RVB(Lattice &lattice,int seed):lattice_(lattice){
    Init(seed); 
  }
 
  // Return the lattice object
  inline Lattice GetLattice(){return lattice_;}

  // Initialization
  void Init(int seed){
    L_        = lattice_.LinSize();
    numSpins_ = lattice_.Nsites();
    ratio_    = 0;
    spins_.resize(numSpins_);
    VB_.resize(numSpins_);
    rgen_.seed(seed);

    // Reset the configuration into the trivial topological sector
    for (int i=0;i<numSpins_;i+=2){
      VB_[i]  = lattice_.neighbours_[i][0];
      VB_[lattice_.neighbours_[i][0]]= i;
    }
    std::uniform_real_distribution<double> distribution(0,1);
    int spin_val;
    // Initialize randomly the spins, according to the dimers
    for (int i=0;i<numSpins_;i+=2){
      spin_val = (distribution(rgen_) < 0.5);
      spins_[i] = spin_val;
      spins_[VB_[i]] = spin_val ^ 1;
    }
  }
  
  // Copy the VB state into the ancillary sytem
  void CopyVBToAncilla(){
    for(int i=0;i<VB_.size()/2;i++){
      VB_[i+VB_.size()/2] = VB_[i]+VB_.size()/2;
    }
  }

  // Set the spin state
  void SetSpins(std::vector<int> spins){
    for(int i=0;i<spins.size();i++){
      spins_[i] = spins[i];
    }
  }

  // Set the valence bond state
  void SetVB(std::vector<int> dimers){
    for(int i=0;i<VB_.size();i++){
      VB_[i] = VB_[i];
    }
  }

  // Set the full state (VB+spins)
  void SetState(RVB &rvb){
    SetSpins(rvb.spins_);
    SetVB(rvb.VB_);
  }

  // Set the regionX for the modified ratio sampling
  void SetRegion(std::vector<int> &region){
    ratio_ = 1;
    regionX_.resize(region.size());
    for(int i=0;i<region.size();i++){
      regionX_[i] = region[i];
    }
  }

  // Perform a local update of VB on a plaquette
  void LocalBondUpdate(int plaq){
    //Check if there are two horizontal VBs in the plaquette
    if(VB_[lattice_.SitesOnPlaquettes_[plaq][0]] == lattice_.SitesOnPlaquettes_[plaq][1] &&
       VB_[lattice_.SitesOnPlaquettes_[plaq][2]] == lattice_.SitesOnPlaquettes_[plaq][3]){
      //Check if spins configurations are compatible
      if(spins_[lattice_.SitesOnPlaquettes_[plaq][0]] != spins_[lattice_.SitesOnPlaquettes_[plaq][3]] &&
         spins_[lattice_.SitesOnPlaquettes_[plaq][1]] != spins_[lattice_.SitesOnPlaquettes_[plaq][2]]) {
        // Perform the update
        VB_[lattice_.SitesOnPlaquettes_[plaq][0]] = lattice_.SitesOnPlaquettes_[plaq][3];
        VB_[lattice_.SitesOnPlaquettes_[plaq][1]] = lattice_.SitesOnPlaquettes_[plaq][2];
        VB_[lattice_.SitesOnPlaquettes_[plaq][2]] = lattice_.SitesOnPlaquettes_[plaq][1];
        VB_[lattice_.SitesOnPlaquettes_[plaq][3]] = lattice_.SitesOnPlaquettes_[plaq][0];
      }
    }
    //Check if there are two vertical VBs in the plaquette
    else if(VB_[lattice_.SitesOnPlaquettes_[plaq][0]] == lattice_.SitesOnPlaquettes_[plaq][3] &&
            VB_[lattice_.SitesOnPlaquettes_[plaq][1]] == lattice_.SitesOnPlaquettes_[plaq][2]){
      //Check if spins configurations are compatible
      if(spins_[lattice_.SitesOnPlaquettes_[plaq][0]] != spins_[lattice_.SitesOnPlaquettes_[plaq][1]] &&
         spins_[lattice_.SitesOnPlaquettes_[plaq][2]] != spins_[lattice_.SitesOnPlaquettes_[plaq][3]]) {
        // Perform the update
        VB_[lattice_.SitesOnPlaquettes_[plaq][0]] = lattice_.SitesOnPlaquettes_[plaq][1];
        VB_[lattice_.SitesOnPlaquettes_[plaq][1]] = lattice_.SitesOnPlaquettes_[plaq][0];
        VB_[lattice_.SitesOnPlaquettes_[plaq][2]] = lattice_.SitesOnPlaquettes_[plaq][3];
        VB_[lattice_.SitesOnPlaquettes_[plaq][3]] = lattice_.SitesOnPlaquettes_[plaq][2];
      }
    }
  }

  // Perform the spin update (given a second compativel copy of the RVB state)
  void SpinUpdate(RVB &rvb){
    std::vector<int> already_updated;     // 1=updated   0=to update
    already_updated.assign(numSpins_,0);
    int next_site,spin_val;
    int num_loops = 0;
    std::uniform_int_distribution<int> distribution(0,1); 
    // If using the ratio sampling, update the spin state such that
    // the overlap <alpha|Swap|beta> is different from zero
    if (ratio_) {
      Swap(regionX_);  // Swap the states
    }
    for(int i=0;i<numSpins_;i++){
      // If spin i has not been updated
      if (already_updated[i] == 0) {
        // Randomly sample a new spin value
        spin_val = distribution(rgen_);
        spins_[i] = spin_val;
        // Spin i hasn now been updated
        already_updated[i] = 1;
        // Go to the next site: this is the site that is in a valence bond state
        // with spin i in the RVB state of this object
        next_site = VB_[i];
        // Keep going through the loop until you find a site that is already 
        // updated, i.e. the initial site
        while(already_updated[next_site] == 0){
          // If the next site is already updated, there must be an error
          if(already_updated[next_site] == 1) {
            std::cout<<"Error - The loop was already updated!"<<std::endl;
            exit(0);
          }
          // Otherwise, flag the next site as updated
          else{ 
            already_updated[next_site] = 1;
          }
          // Flip the spin to keep the correct staggered pattern
          spin_val = spin_val ^ 1;
          spins_[next_site] = spin_val;
          // Go to the next site: this is the site that is in a valence bond state
          // with spin i in the RVB state complementary to the current object
          next_site = rvb.VB_[next_site];
          // If the next site has not been yet updated, do it
          if(already_updated[next_site] != 1)
            already_updated[next_site] = 1;
          // Otherwise, we reached the beginning of the loop, then exit 
          else break;
          // Flip the spin to keep the correct staggered pattern
          spin_val = spin_val ^ 1;
          spins_[next_site] = spin_val;
          // Go to the next site: this is the site that is in a valence bond state
          // // with spin i in the RVB state complementary to the current object
          next_site = VB_[next_site];
        }
        //Count the total number of loops
        num_loops++;
      }
    }
    // Set the spin state of rvb to the one of the current object
    rvb.SetSpins(spins_);
    // Record the number of loops for both RVB states
    num_loops_ = num_loops;
    rvb.num_loops_ = num_loops;
    // If ratio is on, the unswap the states
    if (ratio_) Swap(regionX_);  //unswap the states
  }

  // Calculate the overlap between the object and a different RVB state
  // This is equivalente to the spin update sequence, without sampling
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
            std::cout<<"Error - The loop was already counted!"<<std::endl;
          else 
            already_counted[next_site] = 1;
          next_site = rvb.VB_[next_site];
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

  // Swap the configurations between real and ancillary system
  // according to regionA
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


};

#endif
