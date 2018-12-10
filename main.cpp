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
  int nregions;
  MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
  MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);

  //Parameters
  Parameters pars(totalnodes,mynode);
  pars.ReadParameters(argc,argv);
  pars.geometry_ = "square"; 
  typedef SquareLattice Lattice;
  Lattice lattice(pars.L_);
  lattice.BuildRegions(pars.geometry_,pars.reg_inc_);
  nregions = lattice.Nregions();
  WriteSimulationInfos(pars);
  WriteEntanglementRegions(pars,lattice);
  
  for(int k=0;k<nregions;k++){
     
    if(k==0) pars.ratio_=0;
    else     pars.ratio_=1;
    
    //Create simulation folder
    //SetupSimulationDir(pars);

    //Quantum Monte Carlo
    QMC<Lattice> qmc(lattice,pars);
    qmc.LoadRegion(pars,k);
    qmc.QMCrun();
    
    //Data Analysis
    Stats stats(pars.nMC_);
    stats.SimpleStat(qmc.RenyiEntropy_); 
 
    //Record
    std::string fname = MeasurementName(pars,"renyi");
    //std::cout<<fname<<std::endl;
    std::ofstream fout(fname, std::ios_base::app | std::ios_base::out);
    fout<<k<<" \t";
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
