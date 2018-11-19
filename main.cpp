#include <iostream>

//#include "hypercube.cpp"
#include "qmc.hpp"
#include "stats.hpp"
#include "mpi.h"
#include <ctime>
#include <boost/format.hpp>
//#include "parameters.hpp"

int main(int argc, char* argv[]){
  
  clock_t time_a = std::clock();
  MPI_Init(&argc,&argv);
  int mynode;
  MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
  int totalnodes;
  MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);

  

  Parameters pars(totalnodes,mynode);
  int nsamples_node = (std::ceil(double(pars.nMC_) / double(totalnodes)));
  std::string foutName;
  
  typedef SquareLattice Lattice;
  SquareLattice lattice(pars.L_);
  //lattice.PrintLattice();
  
  QMC<Lattice> qmc(lattice,pars);;    
  qmc.QMCrun();
  

  Stats stats(pars.nMC_);
  
  // Spin-Spin Correlation Function
  foutName = "data/2dSquareRVB_SpinSpinCorrelation_L"+boost::str(boost::format("%d") % pars.L_)+".txt";
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
