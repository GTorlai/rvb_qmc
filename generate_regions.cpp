#include <iostream>

#include <string>
#include <fstream>
#include <vector>
//#include "mpi.h"
//#include <ctime>
//#include <boost/format.hpp>
//#include "parameters.hpp"
#include "square_lattice.hpp"
//#include "cubic_lattice.hpp"
//#include "rvb.hpp"
//#include "qmc.hpp"

int main(){
  
  int L = 16;
  std::string fname;
  
  typedef SquareLattice Lattice;
  //typedef CubicLattice Lattice;
  
  Lattice lattice(L);
  // Cylinder Regions
  for(int w=1;w<L;w++){
    lattice.BuildRegionCylinder(w);
    fname = "regions/regionA_" + std::to_string(w) + ".txt";
    std::ofstream fout(fname);
    lattice.SaveRegion(fout);
    fout.close();
  }
  for(int w=2;w<L;w++){
    lattice.BuildRegionCylinder(w-1);
    fname = "regions/regionX_" + std::to_string(w) + ".txt";
    std::ofstream fout(fname);
    lattice.SaveRegion(fout);
    fout.close();
  }
}
