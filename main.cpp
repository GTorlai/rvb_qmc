#include <iostream>

//#include "hypercube.cpp"
#include "qmc.hpp"
#include "stats.hpp"
#include "mpi.h"
#include <ctime>
#include <boost/format.hpp>

int main(int argc, char* argv[]){
  
  clock_t time_a = std::clock();
  MPI_Init(&argc,&argv);
  int mynode;
  MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
  int totalnodes;
  MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);


  int D = 2;
  int L = 16;

  int seedQMC= 13220+13*mynode;
  int seedBRA= 16382+15*mynode;
  int seedKET= 18209+17*mynode;
  int nburn = 1000;
  int nsamples = 1000000; 
  int nsamples_node = (std::ceil(double(nsamples) / double(totalnodes)));
  std::string foutName;

  typedef SquareLattice Lattice;
  SquareLattice lattice(L);
  //lattice.PrintLattice();
  
  QMC<Lattice> qmc(lattice,seedQMC,seedBRA,seedKET);    
  qmc.QMCrun(nburn,nsamples_node);
  

  Stats stats(nsamples);
  
  // Spin-Spin Correlation Function
  foutName = "data/2dSquareRVB_SpinSpinCorrelation_L"+boost::str(boost::format("%d") % L)+".txt";
  std::ofstream fout(foutName);
  stats.SimpleStat(qmc.SpinSpinCorrelation_);
  stats.SaveVectorStats(fout);
  fout.close();
  
  
  if(mynode == 0){
    clock_t time_b = std::clock();
    unsigned int total_time_ticks = (unsigned int)(time_b - time_a);
    double time_s = float(total_time_ticks)/double(CLOCKS_PER_SEC);
    std::cout<<std::endl<<"Time elapsed: " << time_s << " seconds." << std::endl;
  }
  MPI_Finalize();
}
