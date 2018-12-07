#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <stdio.h>
#include <stdlib.h>

// Parameter Class
class Parameters{

public:
  
  int D_;
  int L_;
  int Wx_;
  int Wy_;
  int nburn_;
  long long int nMC_;
  int ratio_;
  int regionID_;
  int seed_bra_; 
  int seed_ket_;  
  int seed_qmc_;  
  int totalnodes_;
  int mynode_;
  std::string geometry_;

  Parameters(int &totalnodes,int &mynode) {
    D_ = 2;
    L_ = 8;
    Wx_ = 0;
    Wy_ = 0;
    ratio_=1;
    regionID_=1;
    nburn_ = 100000;
    nMC_ = 100000;
    seed_bra_ = 16382+1235*mynode;
    seed_ket_ = 18209+1127*mynode;
    seed_qmc_ = 13220+1663*mynode;
    totalnodes_ = totalnodes;
    mynode_ = mynode;
    geometry_ = "cylinder";
  }
    
  // Read parameters from the command line
  void ReadParameters(int argc,char** argv){
    std::string flag;
    
    flag = "-L";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) L_=atoi(argv[i+1]);
    }
    flag = "-nMC";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) nMC_=atoll(argv[i+1]);
    }
    flag = "-r";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) ratio_=atoi(argv[i+1]);
    }
    flag = "-reg";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) regionID_=atoi(argv[i+1]);
    }
  }
};
#endif
