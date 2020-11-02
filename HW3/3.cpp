#include <iostream>
#include <random>
#include "./LSE_vs_Newton.cpp"
#include "../../lib/matplotlib-cpp-master/matplotlibcpp.h"

using namespace std;

vector<vector<double>> identity_matrix(int n){
  vector<vector<double>> ret;
  for(int i = 0;i < n ;i++){
    vector<double> row(n,0);
    for(int j = 0;j < n ;j++){
      if(i==j){
        row[j] = 1;
      }      
    }
    ret.push_back(row);
  }
  return ret;
}
vector<vector<double>> phi_matrix(int degree, vector<double> xs){
  vector<vector<double>> ret;
  int data_size = xs.size();
  for(int i = 0;i < data_size ;i++){
    vector<double> row;
    for(int j = 0;j < degree ;j++){
      row.push_back(pow(xs[i],j));
    }
    ret.push_back(row);
  }
  return ret;
}


double gaussian_rand(double mean,double deviation){
  mt19937 rng(std::random_device{}());
  uniform_real_distribution<double> distribution(0.0,1.0);
  double U[2] = {};

  for (int i=0; i<2; ++i) {
    U[i] = distribution(rng);
  }

  ////use Box–Muller transform
  double result = sqrt(-2*log(U[0]))*cos(2*M_PI*U[1]);

  return result*deviation+mean;
}

vector<double> gen_linear_data(int n,int a,vector<double> w){

  mt19937 rng(std::random_device{}());
  uniform_real_distribution<double> x_distribution(-1.0,1.0);
  double x = x_distribution(rng);
  double y = 0;
  for(int i = 0; i < n; i++){
    y += w[i]*pow(x,i);
  }
  y += gaussian_rand(0,a);

  vector<double> ret{x,y};
  return ret;
}

namespace plt = matplotlibcpp;
int main()
{
  int a = 1,b = 1,n = 4;
  vector<double> w{1,2,3,4};
  ////plot ground truth
  double min_x = 2,max_x = 2;
  int point_nums = 100;
  std::vector<double> x(point_nums),y(point_nums),y_p_sigma(point_nums),y_m_sigma(point_nums);
  for(int i=0; i < point_nums; i++) {
    x[i] = -2+i/25.0;
    for(int j = 0; j < w.size(); j++){
      y[i] +=  w[j]*pow(x[i],j);
    }
    y_p_sigma[i] = y[i]+1;
    y_m_sigma[i] = y[i]-1;
  }
  plt::subplot(2,2,1);
  plt::plot(x,y);
  plt::plot(x,y_p_sigma,"r");
  plt::plot(x,y_m_sigma,"r");
  ////plot ground truth

  vector<double> xs;
  vector<double> empty_row;
  vector<vector<double>> ts;
  ts.push_back(empty_row);

  int round = 0;
  int cov_count = 0;
  vector<double> init_mean_row(1,0);
  vector<vector<double>> prev_mean(n,init_mean_row);
  vector<vector<double>> prev_covariance = identity_matrix(n);
  while(true){
    // ////test
    // vector<double> new_data{0,0};
    // if(test_count == 1){
    //   new_data = {-0.64152, 0.19039};
    //   test_count++;
    // }else if(test_count == 2){
    //   new_data = {0.07122, 1.63175};
    //   test_count++;
    // }else if(test_count == 3){
    //   new_data = {-0.19330, 0.24507};
    //   test_count++;
    // }else{
    //   new_data = gen_linear_data(n,a,w);
    // }
    round++;
    vector<double> new_data = gen_linear_data(n,a,w);

    ////perdict y of specific x
    double perdict_mean = 0;
    for (int i = 0; i< n; i++){
      perdict_mean += prev_mean[i][0]*pow(new_data[0],i);
    }
    vector<double> input_x{new_data[0]}; 
    double perdict_variance = 1 + matrix_multiply(matrix_multiply(phi_matrix(n,input_x),prev_covariance),matrix_transpose(phi_matrix(n,input_x)))[0][0];
    
    // cout << "Predictive distribution ~ N(" << perdict_mean;
    // cout << ", " << perdict_variance << ")" << endl;
  
    xs.push_back(new_data[0]);
    ts[0].push_back(new_data[1]);
    vector<vector<double>> S_inv = matrix_add(matrix_multiply_constant(identity_matrix(n),a),matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,xs)),phi_matrix(n,xs)),b));
    vector<vector<double>> S = matrix_inverse_byLU(S_inv);
    prev_covariance = S;

    vector<vector<double>> mean = matrix_multiply_constant(matrix_multiply(matrix_multiply(S,matrix_transpose(phi_matrix(n,xs))),matrix_transpose(ts)),b);

    ////Test for if coverage
    bool all_cov = true;
    for (int i = 0; i < prev_mean.size(); i++){
      int a = prev_mean[i][0]*100, b = mean[i][0]*100;
      if ( a != b){
        all_cov = false;
      }
    }

    if(all_cov){
      cov_count++;
    }else{
      cov_count = 0;
    }

    if(round == 10|| round == 50|| cov_count == 3){
      int grid = 0;
      if(round == 10){
        grid = 3;
      }else if(round == 50){
        grid = 4;
      }else{
        grid = 2;
      }
      
      plt::subplot(2,2,grid);
      double min_x = 2,max_x = 2;
      int point_nums = 100;
      std::vector<double> x(point_nums),y(point_nums),y_p_sigma(point_nums),y_m_sigma(point_nums);
      for(int i=0; i < point_nums; i++) {
        x[i] = -2+i/25.0;
        for(int j = 0; j < w.size(); j++){
          y[i] +=  mean[j][0]*pow(x[i],j);
        }
        vector<double> input_x{x[i]}; 
        double variance_at_x = 1 + matrix_multiply(matrix_multiply(phi_matrix(n,input_x),S),matrix_transpose(phi_matrix(n,input_x)))[0][0];
        y_p_sigma[i] = y[i]+variance_at_x;
        y_m_sigma[i] = y[i]-variance_at_x;
      }
      plt::plot(x,y);
      plt::plot(x,y_p_sigma,"r");
      plt::plot(x,y_m_sigma,"r");
      plt::plot(xs,ts[0],"ro");
      if(cov_count == 3){
        plt::show();
        cout << ts[0].size() << endl;
        cout << "stop here" << endl;
    }
    } 
    prev_mean = mean;
  }
  
  cout << "for dubug" << endl;
  return 0;

}