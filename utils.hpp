#ifndef UTILS_HPP
#define UTILS_HPP
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "parameters.hpp"
#include <boost/format.hpp>
//#include <boost/filesystem.hpp>
#include <iostream>
#include "square_lattice.hpp"

std::string SimulationName(Parameters &pars){
  std::string fname;
  //fname  = "../data/sim_";
  //fname = "/Users/gtorlai/Work/Projects/2018/rvb/data/sim";
  fname = "sim";
  fname += boost::str(boost::format("%d") % pars.D_) + "D_";
  fname += "L" + boost::str(boost::format("%d") % pars.L_);
  fname += "_W" + boost::str(boost::format("%d") % pars.Wx_)+boost::str(boost::format("%d") % pars.Wy_);
  return fname;
}

std::string MeasurementName(Parameters &pars,const char* observable){
  std::string fname;
  std::string obs(observable);
  fname  = SimulationName(pars) + "/measurements/";
  fname += obs + "_rank"+ boost::str(boost::format("%d") %pars.mynode_);
  fname += ".txt";
  return fname;
}
void WriteSimulationInfos(Parameters &pars){

  std::string main_path   = SimulationName(pars);
  if(pars.mynode_==0){
    std::string info_file   = main_path + "/sim_infos.txt";
    //Write simulation parameters
    std::ofstream fout(info_file);
    fout<<"Simulation Parameters"<<std::endl<<std::endl;
    fout<<"Dimension = "<<pars.D_<<std::endl;
    fout<<"Linear size = "<<pars.L_<<std::endl;
    fout<<"Topological Sector = ("<<pars.Wx_<<","<<pars.Wy_<<")"<<std::endl;
    fout<<"Burn-in steps = "<<pars.nburn_<<std::endl;
    fout<<"Monte Carlo steps = "<<pars.nMC_<<std::endl;
    fout<<"Entanglement geometry = "<<pars.geometry_<<std::endl;
    if(pars.ratio_) fout<<"Using ratio"<<std::endl;
    else            fout<<"Not using ratio"<<std::endl;
    fout<<"Increment = " <<pars.reg_inc_<<std::endl;
    fout<<"Total number of processors = "<<pars.totalnodes_<<std::endl;
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //Write seeds
  std::string seed_file   = main_path + "/seeds/rank"+boost::str(boost::format("%d") % pars.mynode_)+".txt";
  std::ofstream fout_seed(seed_file);
  fout_seed<<pars.seed_bra_<<"  "<<pars.seed_ket_<<"  "<<pars.seed_qmc_<<std::endl;
  fout_seed.close();
}

template<typename Lattice> void WriteEntanglementRegions(Parameters &pars,Lattice &lattice){
  if(pars.mynode_==0){
    std::string fname;
    std::string path = SimulationName(pars) + "/regions/";
    for(int i=0;i<lattice.regions_.size();i++){
      //fname = SimulationName(pars) + "/regions/";
      fname = path + "regionA_" + std::to_string(i) + ".txt";
      std::ofstream fout(fname);
      for(int y=0;y<pars.L_;y++) {
        for(int x=0;x<pars.L_;x++) {
          fout << lattice.regions_[i][lattice.Index(x,y)] << " ";
        }
        fout<<std::endl;
      }
      fout.close();
    }
    if(pars.ratio_){
      for(int i=0;i<lattice.regions_.size()-1;i++){
        //fname = SimulationName(pars) + "/regions/";
        fname = path + "regionX_" + std::to_string(i+1) + ".txt";
        std::ofstream fout(fname);
        for(int y=0;y<pars.L_;y++) {
          for(int x=0;x<pars.L_;x++) {
            fout << lattice.regions_[i][lattice.Index(x,y)] << " ";
          }
          fout<<std::endl;
        }
        fout.close();
      }
    }
    fname = path + "index_regions.txt";
    std::ofstream fout(fname);
    for(int i=0;i<lattice.regionIndices_.size();i++){
      fout << lattice.regionIndices_[i]<<std::endl;
    }
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


//template<typename Lattice> void WriteEntanglementRegions(Parameters &pars,Lattice &lattice){
//  if(pars.mynode_==0){
//    std::string fname;
//    std::string path = SimulationName(pars) + "/regions/";
//    for(int w=1;w<pars.L_;w++){
//      if(pars.geometry_ == "cylinder")      lattice.BuildRegionCylinder(w);
//      else if (pars.geometry_ == "square")  lattice.BuildRegionRectangle(w,w,0,0);
//      else {
//        std::cout<<"Entanglement region not recognized"<<std::endl;
//        exit(0);
//      }
//      fname = path + "regionA_" + std::to_string(w) + ".txt";
//      std::ofstream fout(fname);
//      for(int y=0;y<pars.L_;y++) {
//        for(int x=0;x<pars.L_;x++) {
//          fout << lattice.regionA_[lattice.Index(x,y)] << " ";
//        }
//        fout<<std::endl;
//      }
//      fout.close();
//    }
//    if(pars.ratio_){
//      for(int w=2;w<pars.L_;w++){
//        if(pars.geometry_ == "cylinder")      lattice.BuildRegionCylinder(w-1);
//        else if (pars.geometry_ == "square")  lattice.BuildRegionRectangle(w-1,w-1,0,0);
//        else {
//          std::cout<<"Entanglement region not recognized"<<std::endl;
//          exit(0);
//        }
//        fname = path + "regionX_" + std::to_string(w) + ".txt";
//        std::ofstream fout(fname);
//        for(int y=0;y<pars.L_;y++) {
//          for(int x=0;x<pars.L_;x++) {
//            fout << lattice.regionA_[lattice.Index(x,y)] << " ";
//          }
//          fout<<std::endl;
//        }
//        fout.close();
//      }
//    }
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//}

//void MakeDir(std::string &dir_name){
//  if ( boost::filesystem::exists(dir_name) ){
//    //std::cout<<"The folder already exists!"<<std::endl;
//    //exit(0);
//  }
//  else {
//    boost::filesystem::create_directory(dir_name);
//  }
//}

//void SetupSimulationDir(Parameters &pars){
//
//  if(pars.mynode_==0){
//    std::string main_path   = SimulationName(pars);
//    std::string regions_dir = main_path + "/regions/";
//    std::string meas_dir    = main_path + "/measurements/";
//    std::string seed_dir    = main_path + "/seeds/";
//    MakeDir(main_path);
//    MakeDir(regions_dir);
//    MakeDir(meas_dir);
//    MakeDir(seed_dir);
//  }
//}



#endif
