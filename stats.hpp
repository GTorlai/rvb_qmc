#ifndef STATS_HPP
#define STATS_HPP

//rvb.cpp
//Class that build the rvb state as a dimer + spin

#include <vector>
#include <iostream>
#include <random>
#include <math.h> 


class Stats{

private:
  
  int Nmeasurements_;
  int Nbins_;

  double scalar_avg_;
  double scalar_var_;
  double scalar_err_;

  std::vector<double> vector_avg_;
  std::vector<double> vector_var_;
  std::vector<double> vector_err_;
  
public:
  //Functions
  Stats(int Nmeasurements):Nmeasurements_(Nmeasurements){
  }

  void Reset(){

  }


  void SimpleStat(std::vector<double> &data){
    scalar_avg_ = 0.0;
    scalar_var_ = 0.0;
    scalar_err_ = 0.0;
    double delta  = 0.0;
    double delta2 = 0.0;
    double m2 = 0.0;

    for(int i=0;i<data.size();i++){
      delta = data[i] - scalar_avg_;
      scalar_avg_ += delta / double(i+1);
      delta2 = data[i] - scalar_avg_;
      m2 += delta * delta2;
    }

    scalar_var_ = m2 / double(data.size() - 1);
    scalar_err_ = sqrt(m2 / double(data.size()*(data.size() - 1)));
    
    //printf("Expectation value = %.10f  +-  %.10f\n",scalar_avg_,scalar_err_);
  }

  void SimpleStat(std::vector<std::vector<double> > &data){
   
    //std::ofstream fout("prova.txt");
    int num_obs = data[0].size();
    vector_avg_.resize(data[0].size());
    vector_var_.resize(data[0].size());
    vector_err_.resize(data[0].size());
    
    std::fill(vector_avg_.begin(),vector_avg_.end(),0.0);

    std::vector<double> delta(data[0].size());
    std::vector<double> delta2(data[0].size());
    std::vector<double> m2(data[0].size());
    std::fill(delta.begin(),delta.end(),0.0);
    std::fill(delta2.begin(),delta2.end(),0.0);
    std::fill(m2.begin(),m2.end(),0.0);

    for(int i=0;i<data.size();i++){
      for(int j=0;j<data[i].size();j++){
        //fout<<data[i][j]<<" ";
        delta[j] = data[i][j] - vector_avg_[j];
        vector_avg_[j] += delta[j] / double(i+1);
        delta2[j] = data[i][j] - vector_avg_[j];
        m2[j] += delta[j] * delta2[j];
      }
      //fout<<std::endl;
    }

    for(int j=0;j<vector_avg_.size();j++){
      vector_var_[j] = m2[j] / double(data.size() - 1);
      vector_err_[j] = sqrt(m2[j] / double(data.size()*(data.size() - 1)));
      printf("Expectation value = %.10f  +-  %.10f\n",vector_avg_[j],vector_err_[j]);
    }
  }
};

#endif
