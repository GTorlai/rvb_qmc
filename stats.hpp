#ifndef STATS_HPP
#define STATS_HPP

#include <vector>
#include <iostream>
#include <random>
#include <math.h> 
#include <mpi.h>
#include <fstream>
class Stats{

private:
  
  int Nmeasurements_;
  int Nbins_;
  int totalnodes_;
  int mynode_;
  double scalar_avg_;
  double scalar_var_;
  double scalar_err_;

  std::vector<double> vector_local_avg_;
  std::vector<double> vector_local_var_;
  std::vector<double> vector_local_err_;
  std::vector<double> vector_avg_;
  std::vector<double> vector_err_;
  
public:
  //Functions
  Stats(int Nmeasurements):Nmeasurements_(Nmeasurements){
    mynode_=0;
    totalnodes_=1;
    //MPI_Comm_size(MPI_COMM_WORLD, &totalnodes_);
    //MPI_Comm_rank(MPI_COMM_WORLD, &mynode_);
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
   
    int num_obs = data[0].size();
    vector_local_avg_.resize(data[0].size());
    vector_local_var_.resize(data[0].size());
    vector_local_err_.resize(data[0].size());
    vector_avg_.resize(data[0].size());
    vector_err_.resize(data[0].size());
    
    std::fill(vector_local_avg_.begin(),vector_local_avg_.end(),0.0);
    std::fill(vector_avg_.begin(),vector_avg_.end(),0.0);
    std::fill(vector_err_.begin(),vector_err_.end(),0.0);

    std::vector<double> delta(data[0].size());
    std::vector<double> delta2(data[0].size());
    std::vector<double> m2(data[0].size());
    std::fill(delta.begin(),delta.end(),0.0);
    std::fill(delta2.begin(),delta2.end(),0.0);
    std::fill(m2.begin(),m2.end(),0.0);

    for(int i=0;i<data.size();i++){
      for(int j=0;j<data[i].size();j++){
        //fout<<data[i][j]<<" ";
        delta[j] = data[i][j] - vector_local_avg_[j];
        vector_local_avg_[j] += delta[j] / double(i+1);
        delta2[j] = data[i][j] - vector_local_avg_[j];
        m2[j] += delta[j] * delta2[j];
      }
      //fout<<std::endl;
    }
    for(int j=0;j<vector_local_avg_.size();j++){
      vector_local_var_[j] = m2[j] / double(data.size() - 1);
      //vector_local_err_[j] = std::sqrt(m2[j] / double(data.size()*(data.size() - 1)));
      vector_local_err_[j] = m2[j] / double(data.size()*(data.size() - 1));
    }
  
    //MPI_Reduce(&vector_local_avg_[0],&vector_avg_[0],vector_avg_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    //MPI_Reduce(&vector_local_err_[0],&vector_err_[0],vector_err_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
    
    //for(int j=0;j<vector_avg_.size();j++){ 
    //  if (mynode_ == 0){
    //    printf("Expectation value = %.10f  +-  %.10f\n",vector_avg_[j]/double(totalnodes_),std::sqrt(vector_err_[j]));
    //  }
    //}
    
    for(int j=0;j<vector_local_avg_.size();j++){ 
      printf("Expectation value = %.10f  +-  %.10f\n",vector_local_avg_[j]/double(totalnodes_),std::sqrt(vector_local_err_[j]));
    }
  }

  //void BinnedStat(std::vector<std::vector<double> > &data){
  // 
  //  int num_obs = data[0].size();
  //  vector_local_avg_.resize(data[0].size());
  //  vector_local_var_.resize(data[0].size());
  //  vector_local_err_.resize(data[0].size());
  //  vector_avg_.resize(data[0].size());
  //  vector_err_.resize(data[0].size());
  //  
  //  std::fill(vector_local_avg_.begin(),vector_local_avg_.end(),0.0);
  //  std::fill(vector_avg_.begin(),vector_avg_.end(),0.0);
  //  std::fill(vector_err_.begin(),vector_err_.end(),0.0);

  //  std::vector<double> delta(data[0].size());
  //  std::vector<double> delta2(data[0].size());
  //  std::vector<double> m2(data[0].size());
  //  std::fill(delta.begin(),delta.end(),0.0);
  //  std::fill(delta2.begin(),delta2.end(),0.0);
  //  std::fill(m2.begin(),m2.end(),0.0);

  //  int num_levels;
  //  int min_bins = 20;
  //  double Delta0,Delta,tau;
  //  num_levels = int(std::floor(std::log2(data.size()/min_bins)));
  //  std::vector<std::vector<double> > bins;
  //  std::vector<std::vector<double> > bins_tmp;
  //  std::vector<std::vector<double> > bin_err;
  //  bins.resize(data.size(),std::vector<double>(data[0].size()));
  //  bins = data;
  //  for(int i=0;i<bins.size();i++){
  //    for(int j=0;j<bins[i].size();j++){
  //      //fout<<data[i][j]<<" ";
  //      delta[j] = bins[i][j] - vector_local_avg_[j];
  //      vector_local_avg_[j] += delta[j] / double(i+1);
  //      delta2[j] = bins[i][j] - vector_local_avg_[j];
  //      m2[j] += delta[j] * delta2[j];
  //    }
  //  }
  //  for(int j=0;j<vector_local_avg_.size();j++){
  //    vector_local_err_[j] = sqrt(m2[j] / double(bins.size()*(bins.size() - 1)));
  //  }
  //  MPI_Reduce(&vector_local_avg_[0],&vector_avg_[0],vector_avg_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //  MPI_Reduce(&vector_local_err_[0],&vector_err_[0],vector_err_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //  
  //  if(mynode_ == 0){
  //    std::cout<<"Bin level = "<< 0 << "   Error = "<<vector_err_[0]/double(totalnodes_)<<std::endl;
  //  } 

  //  for (int l=1;l<num_levels;l++){
  //    if(bins.size() % 2 ==0) {
  //      bins_tmp.resize(int(bins.size()/2),std::vector<double>(bins[0].size()));
  //      for(int i=0;i<bins_tmp.size();i++){
  //        for(int j=0;j<bins_tmp[0].size();j++){
  //          bins_tmp[i][j] = (bins[2*i][j]+bins[2*i+1][j])/2.0;
  //        }
  //      }
  //    }
  //    else{
  //      bins_tmp.resize(int(bins.size()-1)/2),std::vector<double>(bins[0].size());
  //      for(int j=0;j<bins_tmp[0].size();j++){
  //        bins_tmp[0][j] =(bins[0][j]+bins[1][j]+bins[2][j])/3.0;
  //      }
  //      for(int i=1;i<bins_tmp.size();i++){
  //        for(int j=0;j<bins_tmp[0].size();j++){
  //          bins_tmp[i][j] = (bins[2*i+1][j]+bins[2*i+2][j])/2.0;
  //        }
  //      }
  //    }
  //    bins.resize(bins_tmp.size(),std::vector<double>(bins_tmp[0].size()));
  //    bins = bins_tmp;
  //    std::fill(vector_local_avg_.begin(),vector_local_avg_.end(),0.0);
  //    std::fill(vector_avg_.begin(),vector_avg_.end(),0.0);
  //    std::fill(vector_err_.begin(),vector_err_.end(),0.0);
  //    std::fill(delta.begin(),delta.end(),0.0);
  //    std::fill(delta2.begin(),delta2.end(),0.0);
  //    std::fill(m2.begin(),m2.end(),0.0);
  //    
  //    for(int i=0;i<bins.size();i++){
  //      for(int j=0;j<bins[i].size();j++){
  //        //fout<<data[i][j]<<" ";
  //        delta[j] = bins[i][j] - vector_local_avg_[j];
  //        vector_local_avg_[j] += delta[j] / double(i+1);
  //        delta2[j] = bins[i][j] - vector_local_avg_[j];
  //        m2[j] += delta[j] * delta2[j];
  //      }
  //    }
  //    for(int j=0;j<vector_local_avg_.size();j++){
  //      vector_local_err_[j] = sqrt(m2[j] / double(bins.size()*(bins.size() - 1)));
  //    }
  //    MPI_Reduce(&vector_local_avg_[0],&vector_avg_[0],vector_avg_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //    MPI_Reduce(&vector_local_err_[0],&vector_err_[0],vector_err_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //    //MPI_Allreduce(&vector_local_avg_[0],&vector_avg_[0],vector_avg_.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //    //MPI_Allreduce(&vector_local_err_[0],&vector_err_[0],vector_err_.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //    if(mynode_ == 0){
  //      std::cout<<"Bin level = "<< l << "   Error = "<<vector_err_[0]/double(totalnodes_)<<std::endl;
  //    } 
  //  }
  //  //for(int j=0;j<vector_local_avg_.size();j++){
  //  //  vector_local_var_[j] = m2[j] / double(data.size() - 1);
  //  //  vector_local_err_[j] = sqrt(m2[j] / double(data.size()*(data.size() - 1)));
  //  //}
  //
  //  //MPI_Reduce(&vector_local_avg_[0],&vector_avg_[0],vector_avg_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //  //MPI_Reduce(&vector_local_err_[0],&vector_err_[0],vector_err_.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 

  //  //for(int j=0;j<vector_avg_.size();j++){ 
  //  //  if (mynode_ == 0){
  //  //    printf("Expectation value = %.10f  +-  %.10f\n",vector_avg_[j]/double(totalnodes_),vector_err_[j]/double(totalnodes_));
  //  //  }
  //  //}
  //}

  void SaveVectorStats(std::ofstream &fout){
    if (mynode_ == 0){ 
      for(int j=0;j<vector_avg_.size();j++){
        fout<<std::scientific<<vector_avg_[j]/double(totalnodes_)<<"  \t";
        fout<<std::scientific<<std::sqrt(vector_err_[j])<<std::endl;
      }
    }
  }
};

#endif
