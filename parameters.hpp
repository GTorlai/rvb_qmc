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
  int nMC_;
  int ratio_;
  int seed_bra_; 
  int seed_ket_;  
  int seed_qmc_;  
  int totalnodes_;
  
//  Parameters(int &totalnodes,int &mynode) {
  Parameters(){ 
    D_ = 2;
    L_ = 4;
    Wx_ = 0;
    Wy_ = 0;
    ratio_=0;
    nburn_ = 10000;
    nMC_ = 1000000;
    seed_bra_ = 16382;//+15*mynode;
    seed_ket_ = 18209;//+17*mynode;
    seed_qmc_ = 13220;//+13*mynode;
    //totalnodes_ = totalnodes;
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
      if(flag==argv[i]) nMC_=atoi(argv[i+1]);
    }
  }
};
#endif
