#include <iostream>

#include "stats.hpp"
#include "mpi.h"
#include <ctime>
#include <boost/format.hpp>
#include "parameters.hpp"
#include "square_lattice.hpp"
#include "cubic_lattice.hpp"
#include "rvb.hpp"
#include "qmc.hpp"

int main(int argc, char* argv[]){
  
  MPI_Init(&argc,&argv);
  int mynode;
  MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
  int totalnodes;
  MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);
  if(mynode==0){
    std::cout<<std::endl;
  }
  
  Parameters pars(totalnodes,mynode);
  
  pars.ReadParameters(argc,argv); 
  
  typedef SquareLattice Lattice;
  //typedef CubicLattice Lattice;
  
  Lattice lattice(pars.L_);
  //lattice.Print();
  //if(mynode==0){
  //lattice.Print();
  //lattice.PrintRegion();
  //}
  QMC<Lattice> qmc(lattice,pars);;    
  qmc.LoadRegion(pars.regionID_);
  qmc.QMCrun();
  
  Stats stats(pars.nMC_);
  
  //std::cout<<" Spin-Spin Correlation Function " << std::endl<<std::endl;
  //stats.SimpleStat(qmc.SpinSpinCorrelation_);
  //for(int j=0;j<stats.vector_local_avg_.size();j++){ 
  //  printf("Expectation value = %.10f  +-  %.10f\n",stats.vector_local_avg_[j]/double(stats.totalnodes_),std::sqrt(stats.vector_local_err_[j]));
  //}
  //std::cout<<std::endl;
  //std::cout<<" Renyi Entanglement Entropy " << std::endl<<std::endl;
  //stats.SimpleStat(qmc.RenyiEntropy_);
  //for(int j=0;j<stats.vector_local_avg_.size();j++){ 
  //  printf("Expectation value = %.10f  +-  %.10f\n",(stats.vector_local_avg_[j]/double(stats.totalnodes_)),std::sqrt(stats.vector_local_err_[j]));
  //}
  //std::cout<<std::endl;

  stats.SimpleStat(qmc.RenyiEntropy_);
  //printf("Expectation value = %.10f  +-  %.10f\n",(stats.scalar_avg_/double(stats.totalnodes_)),stats.scalar_err_);
  std::string fname;
  fname = "renyi"+std::to_string(pars.regionID_)+".txt";
  std::ofstream fout(fname);
  stats.SaveScalarStats(fout);
  //fout.close();
  //}
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
  MPI_Finalize();

}
