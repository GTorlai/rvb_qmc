#include <iostream>

#include "stats.hpp"
#include "mpi.h"
#include <ctime>
#include <boost/format.hpp>
#include "parameters.hpp"
#include "square_lattice.hpp"
#include "rvb.hpp"
#include "qmc.hpp"

int main(int argc, char* argv[]){
  
  Parameters pars;
  pars.ReadParameters(argc,argv); 
  
  typedef SquareLattice Lattice;
  SquareLattice lattice(pars.L_);
  //lattice.Print();
  //if(mynode==0){
  //  lattice.PrintLattice();
  //}
  QMC<Lattice> qmc(lattice,pars);;    
  qmc.QMCrun();
  
  Stats stats(pars.nMC_);
  stats.SimpleStat(qmc.SpinSpinCorrelation_);


  //// Spin-Spin Correlation Function
  //foutName = "data/2dSquareRVB_SpinSpinCorrelation_L"+boost::str(boost::format("%d") % pars.L_);
  //foutName += "_nMC"+boost::str(boost::format("%e") % pars.nMC_) +".txt";
  //std::ofstream fout(foutName);
  //stats.SimpleStat(qmc.SpinSpinCorrelation_);
  ////stats.BinnedStat(qmc.SpinSpinCorrelation_);
  //stats.SaveVectorStats(fout);
  //fout.close();
  //
  //if(mynode == 0){
  //  clock_t time_b = std::clock();
  //  unsigned int total_time_ticks = (unsigned int)(time_b - time_a);
  //  double time_s = float(total_time_ticks)/double(CLOCKS_PER_SEC);
  //  std::cout<<std::endl<<"Time elapsed: " << time_s << " seconds." << std::endl;
  //}
  //MPI_Finalize();

}
