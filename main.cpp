#include <iostream>

//#include "hypercube.cpp"
#include "qmc.hpp"

int main(int argc, char* argv[]){

    int D = 2;
    int L = 4;

    
    typedef SquareLattice Lattice;
    SquareLattice lattice(L);
    lattice.PrintLattice();
    //Hypercube lattice(L,D);  
    //lattice.Print();
    //
    //RVB<Lattice> rvb(lattice);
    //rvb.PrintSpins();
    
    int seedBRA= 16382;
    int seedKET= 18209;
    QMC<Lattice> qmc(lattice,seedBRA,seedKET);    
    qmc.Sweep(10000);
    //qmc.Test(1); 
    
    //////---- PARAMETERS ----//
    //qst::Parameters par;    //Set default initial parameters

    ////Read simulation parameters from command line
    //par.ReadParameters(argc,argv);    //Read parameters from the command line
    ////par.PrintParameters();            //Print parameter on screen
   
    ////---- SPECIFIC PARAMETERS ----/
    //
    //// POSISTIVE - TFIM1d with 10 SPINS
    ////typedef qst::WavefunctionPositive NNState;       //Positive Wavefunction
    ////par.basis_ = "std";
    ////par.nv_=10;
    ////par.nh_=10;
    ////std::string model = "tfim1d_N"+boost::str(boost::format("%d") % par.nv_);
    ////std::string baseName = "../data/"+model+"_";
    ////par.ns_= 10000;
    ////par.nc_= 100;
    //////std::string parName = "parameters_train_benchmark_real.txt";
    //
    ////// REAL - HEISENBERG with 10 spins
    ////typedef qst::WavefunctionReal NNState;       //Complex Wavefunction
    ////par.basis_ = "xy2";
    //////std::string model = "2qubits";
    //////std::string baseName = "../../examples/2qubits_complex/"+model+"_";
    //////std::string model = "heis1d_N4";
    //////std::string baseName = "../../examples/heisenberg_model/"+model+"_";
    //////std::string model = "j1j21d_N4";
    //////std::string baseName = "../../examples/j1j2_model/"+model+"_";
    ////std::string baseName = "../../examples/03_data_generation/";
    ////par.nv_=2;
    ////par.nh_=par.nv_;
    ////par.ns_=1000;
    ////par.nc_=100;
    //////std::string parName = "parameters_benchmark_complex.txt";
    //
    //
    ////// COMPLEX 
    //typedef qst::WavefunctionComplex NNState;       //Complex Wavefunction
    //par.basis_ = "xy1";
    //std::string model = "qubits";
    //std::string baseName = "../data/"+model+"_";
    ////std::string model = "heis1d_N2";
    ////std::string baseName = "../../examples/template/"+model+"_";
    ////std::string baseName = "../../examples/tmp/heisenberg_model/"+model+"_";
    /////std::string model = "j1j21d_N4";
    /////std::string baseName = "../../examples/j1j2_model/"+model+"_";
    //par.nv_=2;
    ////par.nh_=par.nv_;
    //par.ns_=100;
    //par.nc_=10;
    ////std::string parName = "parameters_benchmark_complex.txt";
    //
    //typedef qst::Sgd Optimizer;                     //Stochastic gradient descent
    ////typedef qst::AdaDelta Optimizer;                     //Stochastic gradient descent
    ////typedef qst::Ngd Optimizer;                     //Stochastic gradient descent
    //typedef qst::ObserverPSI<NNState> Observer;              //Observer for Wavefunction

    //////Load the data
    //std::string fileName; 
    //qst::SetNumberOfBases(par);
    //Eigen::VectorXcd target_psi;                //Target wavefunction
    //std::vector<Eigen::VectorXcd> rotated_target_psi;       //Vector with the target wavefunctions in different basis
    //std::vector<std::vector<std::string> > basisSet;        //Set of bases available
    //std::map<std::string,Eigen::MatrixXcd> UnitaryRotations;//Container of the of 1-local unitary rotations
    //Eigen::MatrixXd training_samples(par.ns_,par.nv_);      //Training samples matrix
    //std::vector<std::vector<std::string> > training_bases;  //Training bases matrix

    ////Load data
    //qst::GenerateUnitaryRotations(UnitaryRotations);        //Generate the unitary rotations
    //fileName = baseName + "psi.txt";            
    //qst::LoadWavefunction(par,fileName,target_psi,rotated_target_psi);
    //fileName = baseName + "bases.txt";
    //qst::LoadBasesConfigurations(par,fileName,basisSet);                //Load training samples
    //qst::LoadTrainingData(baseName,par,training_samples,training_bases);//Load training bases
    //////---- OPTIMIZER ----//
    //Optimizer opt(par);         //Construc the optimizer object
    //////---- NEURAL NETWORK STATE ----//
    //NNState nn(par);
    //nn.InitRandomPars(12345,par.w_);
    ////nn.LoadWeights(parName);    
    ////nn.PrintParameters();
    ////---- OBSERVER ----//
    //Observer obs(nn,par.basis_);
    //obs.setWavefunction(target_psi);
    //if (par.basis_.compare("std")!=0){
    //    obs.setBasisRotations(UnitaryRotations);
    //    obs.setBasis(basisSet);
    //    obs.setRotatedWavefunctions(rotated_target_psi);
    //} 
    //
    /////---- TOMOGRAPHY ----//
    ////qst::Tomography<NNState,Observer,Optimizer> tomo(opt,nn,obs,par);
    ////tomo.setBasisRotations(UnitaryRotations);
    ////tomo.Run(training_samples,training_bases);
   

    ////---- TEST ----// 
    //par.PrintParameters();            //Print parameter on screen
    //qst::Test<NNState,Observer> test(nn,obs,par);
    //test.setWavefunction(target_psi);
    //if (par.basis_.compare("std")!=0){
    //    test.setBasisRotations(UnitaryRotations);
    //    test.setBasis(basisSet);
    //    test.setRotatedWavefunctions(rotated_target_psi);
    //}
    //test.RunDerCheck(par.nh_,training_samples,training_bases,1e-8,par.nc_);
}
