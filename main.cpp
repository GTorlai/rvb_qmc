#include <iostream>
#include <iomanip>
#include "stats.hpp"
#include "mpi.h"
#include <ctime>
#include <boost/format.hpp>
#include "parameters.hpp"
#include "square_lattice.hpp"
#include "cubic_lattice.hpp"
#include "rvb.hpp"
#include "qmc.hpp"
#include "utils.hpp"
int main(int argc, char* argv[]){
  //clock_t time_a = std::clock(); 
  MPI_Init(&argc,&argv);
  int mynode;
  int totalnodes;
  MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
  MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);

  //Parameters
  Parameters pars(totalnodes,mynode);
  pars.ReadParameters(argc,argv);
  std::cout<<pars.nMC_<<std::endl;
  for(int k=1;k<pars.L_;k++){
    
    if(k==1) pars.ratio_=0;
    else     pars.ratio_=1;
    pars.regionID_ = k;
    //Create simulation folder
    //SetupSimulationDir(pars);
    WriteSimulationInfos(pars);

    //Lattice 
    typedef SquareLattice Lattice;
    //typedef CubicLattice Lattice;
    Lattice lattice(pars.L_);
    WriteEntanglementRegions(pars,lattice); 
    
    //Quantum Monte Carlo
    QMC<Lattice> qmc(lattice,pars);
    qmc.LoadRegion(pars);
    qmc.QMCrun();
    
    //Data Analysis
    Stats stats(pars.nMC_);
    stats.SimpleStat(qmc.RenyiEntropy_); 
 
    //Record
    std::string fname = MeasurementName(pars,"renyi");
    std::ofstream fout(fname, std::ios_base::app | std::ios_base::out);
    fout<<pars.regionID_<<" \t";
    fout<<std::setprecision(10)<<stats.scalar_local_avg_<<" \t";
    fout<<std::setprecision(10)<<stats.scalar_local_var_<<" \t";
    fout<<stats.size_<<std::endl;
    fout.close();
  }
  //if(mynode == 0){
  //  clock_t time_b = std::clock();
  //  unsigned int total_time_ticks = (unsigned int)(time_b - time_a);
  //  double time_s = float(total_time_ticks)/double(CLOCKS_PER_SEC);
  //  std::cout<<std::endl<<"Region # = " << pars.regionID_ << "  -  Time elapsed: " << time_s << " seconds." << std::endl;
  //}
  MPI_Finalize();
}
