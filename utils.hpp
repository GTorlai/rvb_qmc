#ifndef UTILS_HPP
#define UTILS_HPP
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "parameters.hpp"
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include "square_lattice.hpp"

std::string SimulationName(Parameters &pars){
  std::string fname;
  //fname  = "../data/sim_";
  fname = "/Users/gtorlai/Work/Projects/2018/rvb/data/sim";
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


void MakeDir(std::string &dir_name){
  if ( boost::filesystem::exists(dir_name) ){
    //std::cout<<"The folder already exists!"<<std::endl;
    //exit(0);
  }
  else {
    boost::filesystem::create_directory(dir_name);
  }
}

void SetupSimulationDir(Parameters &pars){

  if(pars.mynode_==0){
    std::string main_path   = SimulationName(pars);
    std::string regions_dir = main_path + "/regions/";
    std::string meas_dir    = main_path + "/measurements/";
    std::string seed_dir    = main_path + "/seeds/";
    MakeDir(main_path);
    MakeDir(regions_dir);
    MakeDir(meas_dir);
    MakeDir(seed_dir);
  }
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

//void WriteEntanglementRegions(Parameters &pars,Lattice &lattice){
template<typename Lattice> void WriteEntanglementRegions(Parameters &pars,Lattice &lattice){
  if(pars.mynode_==0){
    std::string fname;
    std::string path = SimulationName(pars) + "/regions/";
    for(int w=1;w<pars.L_;w++){
      if(pars.geometry_ == "cylinder")      lattice.BuildRegionCylinder(w);
      else if (pars.geometry_ == "square")  lattice.BuildRegionRectangle(w,w,0,0);
      else {
        std::cout<<"Entanglement region not recognized"<<std::endl;
        exit(0);
      }
      fname = path + "regionA_" + std::to_string(w) + ".txt";
      std::ofstream fout(fname);
      for(int y=0;y<pars.L_;y++) {
        for(int x=0;x<pars.L_;x++) {
          fout << lattice.regionA_[lattice.Index(x,y)] << " ";
        }
        fout<<std::endl;
      }
      //SaveRegion(fout);
      fout.close();
    }
    if(pars.ratio_){
      for(int w=2;w<pars.L_;w++){
        if(pars.geometry_ == "cylinder")      lattice.BuildRegionCylinder(w-1);
        else if (pars.geometry_ == "square")  lattice.BuildRegionRectangle(w-1,w-1,0,0);
        else {
          std::cout<<"Entanglement region not recognized"<<std::endl;
          exit(0);
        }
        fname = path + "regionX_" + std::to_string(w) + ".txt";
        std::ofstream fout(fname);
        for(int y=0;y<pars.L_;y++) {
          for(int x=0;x<pars.L_;x++) {
            fout << lattice.regionA_[lattice.Index(x,y)] << " ";
          }
          fout<<std::endl;
        }
        //SaveRegion(fout);
        fout.close();
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}



//std::string TrainingDataName(Parameters & p,std::string &source){
//
//    std::string fileName;
//    fileName = "../data/datasets/";
//    fileName += source;
//    //fileName += "/N" + boost::str(boost::format("%d") % p.nsites_);
//    fileName += "/" + p.model_ + "_N";
//    fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "_ind";
//    fileName += boost::str(boost::format("%d") % p.h_); 
//    //if (source == "ed_noise"){
//    //    fileName += "_pR" + boost::str(boost::format("%.2f") % p.pR_);
//    //    if (p.pR_ == 0.1){
//    //        fileName += "_pR" + boost::str(boost::format("%.1f") % p.pR_);
//    //    }
//    //    fileName += "_pG" + boost::str(boost::format("%.1f") % p.pG_);
//    //}
//    //if (source == "amplitude_damping_ed"){
//    //    fileName += "_pR" + boost::str(boost::format("%.2f") % p.pR_);
//    //    fileName += "_pG" + boost::str(boost::format("%.1f") % p.pG_);
//    //}
//    fileName += "_"+p.basis_+"basis";
//    fileName += "_train.txt"; 
//    //std::cout << fileName << std::endl;
//    return fileName;
//}
//
//std::string WavefunctionName(Parameters & p){
//    std::string fileName;
//    fileName = "../data/wavefunctions";//N";
//    //fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "/wavefunction_"+p.model_+"_N";
//    fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "_ind" + boost::str(boost::format("%d") % p.h_);
//    fileName += ".txt";
//    //std::cout<<fileName<<std::endl;
//    return fileName;
//}
//
//std::string NetworkName(Parameters & p,std::string &source){
//    std::string fileName;
//
//    fileName = "rbmState_rydberg_";
//    fileName += source;
//    fileName += "_N" + boost::str(boost::format("%d") % p.nsites_);
//    fileName += "_nh" + boost::str(boost::format("%d") % p.nh_);
//    fileName += "_" + p.alg_;
//    fileName += boost::str(boost::format("%d") % p.cd_); 
//    fileName += "_nc" + boost::str(boost::format("%d") % p.nc_);
//    if (p.alg_ == "PT"){
//        fileName += "_nrep" + boost::str(boost::format("%d") % p.nrep_);
//    }
//    fileName += "_" + p.opt_;
//    if (p.opt_ == "sgd"){
//        fileName += "_lr";
//        if (p.lr_ > 0.09) {
//            fileName += boost::str(boost::format("%.1f") % p.lr_);
//        }
//        else if (p.lr_ > 0.009) {
//            fileName += boost::str(boost::format("%.2f") % p.lr_);
//        }
//        else if (p.lr_ > 0.0009) {
//            fileName += boost::str(boost::format("%.3f") % p.lr_);
//        }
//        else if (p.lr_ > 0.00009) {
//            fileName += boost::str(boost::format("%.4f") % p.lr_);
//        }
//        else if (p.lr_ > 0.000009) {
//            fileName += boost::str(boost::format("%.5f") % p.lr_);
//        }
//        fileName += "_reg";
//        if (p.l2_ > 0.09) {
//            fileName += boost::str(boost::format("%.1f") % p.lr_);
//        }
//        else if (p.l2_ > 0.009) {
//            fileName += boost::str(boost::format("%.2f") % p.l2_);
//        }
//        else if (p.l2_ > 0.0009) {
//            fileName += boost::str(boost::format("%.3f") % p.l2_);
//        }
//        else if (p.l2_ > 0.00009) {
//            fileName += boost::str(boost::format("%.4f") % p.l2_);
//        }
//        else if (p.l2_ > 0.000009) {
//            fileName += boost::str(boost::format("%.5f") % p.l2_);
//        }
//    }
//    else if (p.opt_ == "adaDelta"){
//        fileName += "_eps1e-5";
//    }
//    fileName += "_bs" + boost::str(boost::format("%d") % p.bs_);
//    fileName += "_w";
//    if (p.w_ > 0.09) {
//        fileName += boost::str(boost::format("%.1f") % p.w_);
//    }
//    else if (p.w_ > 0.009) {
//        fileName += boost::str(boost::format("%.2f") % p.w_);
//    }
//    else if (p.w_ > 0.0009) {
//        fileName += boost::str(boost::format("%.3f") % p.w_);
//    }
// 
//    if (p.opt_ == "ngd"){
//        fileName += "_lambda";
//        if (p.lambda_ > 0.009) {
//            fileName += boost::str(boost::format("%.2f") % p.lambda_);
//        }
//        else if (p.lambda_ > 0.0009) {
//            fileName += boost::str(boost::format("%.3f") % p.lambda_);
//        }
//        else if (p.lambda_ > 0.00009) {
//            fileName += boost::str(boost::format("%.4f") % p.lambda_);
//        }
//    }
//    fileName += "_ns";
//    fileName += boost::str(boost::format("%d") % p.ns_);
//    fileName += "_"+p.basis_+"basis";
//    fileName += "_ind" + boost::str(boost::format("%d") % p.h_);
//    //if (source == "ed_noise"){
//    //    fileName += "_pR" + boost::str(boost::format("%.2f") % p.pR_);
//    //    if (p.pR_ == 0.1){
//    //        fileName += "_pR" + boost::str(boost::format("%.1f") % p.pR_);
//    //    }
//    //    fileName += "_pG" + boost::str(boost::format("%.1f") % p.pG_);
//    //}
//    //if (source == "amplitude_damping_ed"){
//    //    fileName += "_pR" + boost::str(boost::format("%.2f") % p.pR_);
//    //    fileName += "_pG" + boost::str(boost::format("%.1f") % p.pG_);
//    //}
//    return fileName;
//}
//
//std::string RbmWeightsName(Parameters & p,std::string &source){
//    std::string fileName;
//    fileName = "../data/weights/"+source;//+"/N";
//    //fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "/" + NetworkName(p,source);
//    fileName += "_weights.txt";
//    return fileName;
//}
//
//std::string RbmWavefunctionName(Parameters & p,std::string &source){
//    std::string fileName;
//    fileName = "../data/rbmWavefunctions/";
//    fileName += source;
//    //fileName += "/N";
//    //fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "/" + NetworkName(p,source);
//    fileName += "_wavefunction.txt";
//    return fileName;
//}
//
//std::string InfosName(Parameters & p,std::string &source){
//    std::string fileName;
//    fileName = "../data/datasets/";
//    fileName += source;
//    //fileName += "/N" + boost::str(boost::format("%d") % p.nsites_);
//    fileName += "/" + p.model_ + "_N";
//    fileName += boost::str(boost::format("%d") % p.nsites_);
//    fileName += "_"+p.basis_+"basis";
//    fileName += "_infos.txt";
//    //std::cout<<fileName<<std::endl;
//    return fileName;
//}
//
//std::string ObservableName(Parameters & p, std::string & type,std::string &source,std::string &obs){
//    std::string fileName;
//    fileName = "../data/observables";
//    //fileName += "/N" + boost::str(boost::format("%d") % p.nsites_);
//    fileName += type + source;
//    fileName += "N" + boost::str(boost::format("%d") % p.nsites_);
//    return fileName;
//}
//std::string RbmObserverName(Parameters & p,std::string &source){
//    std::string fileName;
//    fileName = "../data/observers/";
//    fileName += source + "/";
//    //fileName += "/N"+ boost::str(boost::format("%d") %p.nsites_) +"/";
//    fileName += NetworkName(p,source);
//    fileName += "_observer.txt";
//    return fileName;
//}
//
//
//}

#endif
