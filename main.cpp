#include <iostream>

//#include "hypercube.cpp"
#include "qmc.hpp"
#include "stats.hpp"
int main(int argc, char* argv[]){

    int D = 2;
    int L = 8;

    typedef SquareLattice Lattice;
    SquareLattice lattice(L);
    //lattice.PrintLattice();
    
    int seedQMC= 13220;
    int seedBRA= 16382;
    int seedKET= 18209;
    int nsweeps = 1000; 
    QMC<Lattice> qmc(lattice,seedQMC,seedBRA,seedKET);    
    qmc.QMCrun(nsweeps/100,nsweeps);
  
//    for(int i=0;i<lattice.LinSize()-1;i++){
//      std::cout<< "<S0 S"<<i+1<<"> = ";
//      printf("%.6f  +-  %.6f\n",qmc.SpinSpinCorrelation_Avg_[i],sqrt(qmc.SpinSpinCorrelation_Var_[i]/float(nsweeps*(nsweeps+1))));
//    }
    std::cout<<std::endl<<std::endl;

    Stats stats(nsweeps);
    //stats.SimpleScalarStat(qmc.NNcorrelation_);
    stats.SimpleStat(qmc.SpinSpinCorrelation_);
}
