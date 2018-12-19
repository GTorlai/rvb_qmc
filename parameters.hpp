#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <stdio.h>
#include <stdlib.h>

// Parameter Class
class Parameters{

public:
  
  int D_;                 // Lattice dimension
  int L_;                 // Linea size of the lattice
  int Wx_;                // Topological sector in X
  int Wy_;                // Topological sector in Y
  int nburn_;             // Number of burn-in steps
  long long int nMC_;     // Number of measurements per bin
  int nbins_;             // Number of bins
  int ratio_;             // Ration flag
  int reg_inc_;           // Region increment
  int num_reg_;           // Number of regions
  int seed_bra_;          // Seed of the bra       
  int seed_ket_;          // Seed of the ket
  int seed_qmc_;          // Seed of the qmc
  int totalnodes_;        // Number of processors
  int mynode_;            // Processor ID
  std::string geometry_;  // Entanglement geometry
  int sim_id_; 

  Parameters(int &totalnodes,int &mynode) {
    D_ = 2;
    L_ = 8;
    Wx_ = 0;
    Wy_ = 0;
    ratio_=1;
    reg_inc_ = 1;
    nburn_ = 100000;
    nMC_ = 100000;
    seed_bra_ = 16382+1235*mynode;
    seed_ket_ = 18209+1127*mynode;
    seed_qmc_ = 13220+1663*mynode;
    totalnodes_ = totalnodes;
    mynode_ = mynode;
    nbins_ = totalnodes_;
    sim_id_ = 1;
    geometry_ = "cylinder";
  }
    
  // Read parameters from the command line
  void ReadParameters(int argc,char** argv){
    std::string flag;
    
    flag = "-L";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) L_=atoi(argv[i+1]);
    }
    flag = "-mc";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) nMC_=atoll(argv[i+1]);
    }
    flag = "-bins";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) nbins_=atoi(argv[i+1]);
    }
    flag = "-reg";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) reg_inc_=atoi(argv[i+1]);
    }
    flag = "-id";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) sim_id_=atoi(argv[i+1]);
    }
  }
};
#endif
