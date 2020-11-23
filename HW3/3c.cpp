#include <iostream>
#include <random>
#include "./LSE_vs_Newton.cpp"
#include "../../lib/matplotlib-cpp-master/matplotlibcpp.h"

using namespace std;

void print_matrix(vector<vector<double>> mat){
  for(int i = 0;i< mat.size();i++){
    for(int j = 0;j< mat[0].size();j++){
      cout << mat[i][j] << ",";
    }  
    cout << endl;
  }
}

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
  int a = 1,n = 4;
  double b = 1;
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
    y_p_sigma[i] = y[i]+a;
    y_m_sigma[i] = y[i]-a;
  }
  plt::subplot(2,2,1);
  plt::title("Ground truth");
  plt::plot(x,y,"k");
  plt::plot(x,y_p_sigma,"r");
  plt::plot(x,y_m_sigma,"r");
  ////plot ground truth

  vector<double> xs;
  vector<double> empty_row;
  vector<vector<double>> ts;
  ts.push_back(empty_row);

  double round = 0;
  int cov_count = 0;
  int test_count = 1;
  double mean_guess = 0, varriance = 0,s = 0;
  int mean_cov = 0,deviation_cov = 0;
  vector<double> init_mean_row(1,0);
  vector<vector<double>> prev_mean(n,init_mean_row);
  vector<vector<double>> prev_covariance = matrix_multiply_constant(identity_matrix(n),1/b);
  vector<vector<double>> S0_inv,mean0;
  while(true){
    round++;
    ////use sample data
    vector<double> new_data{0,0};
    if(test_count == 1){
      new_data = {-0.64152, 0.19039};
      test_count++;
    }else if(test_count == 2){
      new_data = {0.07122, 1.63175};
      test_count++;
    }else if(test_count == 3){
      new_data = {-0.19330, 0.24507};
      test_count++;
    }else{
      new_data = gen_linear_data(n,a,w);
    }
    // vector<double> new_data = gen_linear_data(n,a,w);

    ////count mean & variance overall
    double precision = 1;
    if(round >1){
      double xn = new_data[1];
      // for (int i = 0;i< n;i++){
      //   xn +=pow(new_data[1],i);
      // }
      double meann = ( mean_guess * (round-1) + xn )/( round );
      s = s + ( xn - mean_guess )*( xn - meann );
      mean_guess = meann;
      varriance = s/round;
      precision = varriance;
      // precision = precision + pow(mean_guess,2) - pow(meann,2) + ((pow(xn,2) - precision - pow(mean_guess,2)) / (round+1));
      // mean_guess = meann;
    }else{
      mean_guess = new_data[1];
    }

    ////perdict y of specific x
    double perdict_mean = 0;
    for (int i = 0; i< n; i++){
      perdict_mean += prev_mean[i][0]*pow(new_data[0],i);
    }
    vector<double> input_x{new_data[0]}; 
    double perdict_variance = 1/1 + matrix_multiply(matrix_multiply(phi_matrix(n,input_x),prev_covariance),matrix_transpose(phi_matrix(n,input_x)))[0][0];
    
  
    xs.push_back(new_data[0]);
    ts[0].push_back(new_data[1]);
    vector<vector<double>> mean;
    vector<vector<double>> S;
    // if(round == 1){
      S0_inv = matrix_add(matrix_multiply_constant(identity_matrix(n),b),matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,xs)),phi_matrix(n,xs)),1));
      S = matrix_inverse_byLU(S0_inv);
      prev_covariance = S;
      mean = matrix_multiply_constant(matrix_multiply(matrix_multiply(S,matrix_transpose(phi_matrix(n,xs))),matrix_transpose(ts)),1);      
      mean0 = mean;
    // }else{
    //   vector<vector<double>> S_inv = matrix_add(S0_inv,matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,xs)),phi_matrix(n,xs)),precision));
    //   S = matrix_inverse_byLU(S_inv);
    //   prev_covariance = S;
    //   mean = matrix_multiply(S,matrix_add(matrix_multiply(S0_inv,mean0),matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,xs)),matrix_transpose(ts)),precision)));     
    //   S0_inv = S_inv;
    //   mean0 = mean;
    // }

    cout << "Add data point (" << new_data[0] << ", " << new_data[1] << "):"<< endl << endl;
    cout << "Postirior mean:" << endl;
    print_matrix(mean);
    cout << endl;
    cout << "Postirior variance:" << endl;
    print_matrix(S);
    cout << endl;
    cout << "Predictive distribution ~ N(" << perdict_mean;
    cout << ", " << perdict_variance << ")" << endl;
    cout << "--------------------------------------------------" << endl;

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
        plt::subplot(2,2,3);
        plt::title("After 10 incomes");
      }else if(round == 50){
        plt::subplot(2,2,4);
        plt::title("After 50 incomes");
      }else{
        plt::subplot(2,2,2);
        plt::title("Perdict result");
      }
      double min_x = 2,max_x = 2;
      int point_nums = 100;
      std::vector<double> x(point_nums),y(point_nums),y_p_sigma(point_nums),y_m_sigma(point_nums);
      for(int i=0; i < point_nums; i++) {
        x[i] = -2+i/25.0;
        for(int j = 0; j < w.size(); j++){
          y[i] +=  mean[j][0]*pow(x[i],j);
        }
        vector<double> input_x{x[i]}; 
        double variance_at_x = 1/precision + matrix_multiply(matrix_multiply(phi_matrix(n,input_x),S),matrix_transpose(phi_matrix(n,input_x)))[0][0];
        y_p_sigma[i] = y[i]+variance_at_x;
        y_m_sigma[i] = y[i]-variance_at_x;
      }
      plt::plot(x,y,"k");
      plt::plot(x,y_p_sigma,"r");
      plt::plot(x,y_m_sigma,"r");
      plt::plot(xs,ts[0],"bo");
      if(cov_count == 3){
        plt::show();
        cout << ts[0].size() << endl;
        cout << "stop here" << endl;
        break;
    }
    } 
    prev_mean = mean;
  }
  
  cout << "for dubug" << endl;
  return 0;

}