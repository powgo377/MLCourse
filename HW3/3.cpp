#include <iostream>
#include <random>
#include "./LSE_vs_Newton.cpp"

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

int main()
{
  int a = 1,b = 1,n = 4;
  vector<double> w{1,2,3,4};

  // vector<double> new_data = gen_linear_data(n,a,w);
  // vector<double> test_x{new_data[0]};
  vector<double> xs;
  vector<double> empty_row;
  vector<vector<double>> ts;
  ts.push_back(empty_row);

  // vector<vector<double>> S_inv = matrix_add(matrix_multiply_constant(identity_matrix(n),a),matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,test_x)),phi_matrix(n,test_x)),b));
  // vector<vector<double>> S = matrix_inverse_byLU(S_inv);
  // vector<vector<double>> mean = matrix_multiply_constant(matrix_multiply(matrix_multiply(S,matrix_transpose(phi_matrix(n,test_x))),matrix_transpose(t)),b);


  int cov_count = 0;
  vector<vector<double>> prev_mean;
  while(true){
    vector<double> new_data = gen_linear_data(n,a,w);
  
    xs.push_back(new_data[0]);
    // vector<double> row{test_2[1]};
    ts[0].push_back(new_data[1]);
    // vector<vector<double>> phi = phi_matrix(4,test_x);
    // vector<vector<double>> bphitphi = matrix_multiply_constant(matrix_multiply(matrix_transpose(phi),phi),b);
    // vector<vector<double>> S_N_inv = matrix_add(S_0_inv,bphitphi);
    // vector<vector<double>> S_N = matrix_inverse_byLU(S_N_inv);
    vector<vector<double>> S_inv = matrix_add(matrix_multiply_constant(identity_matrix(n),a),matrix_multiply_constant(matrix_multiply(matrix_transpose(phi_matrix(n,xs)),phi_matrix(n,xs)),b));
    vector<vector<double>> S = matrix_inverse_byLU(S_inv);

    // vector<vector<double>> phit = matrix_transpose(phi_matrix(n,test_x));
    // vector<vector<double>> sphit = matrix_multiply(S,matrix_transpose(phi_matrix(n,test_x)));
    // vector<vector<double>> sphitt= matrix_multiply(matrix_multiply(S,matrix_transpose(phi_matrix(n,test_x))),matrix_transpose(t));
    vector<vector<double>> mean = matrix_multiply_constant(matrix_multiply(matrix_multiply(S,matrix_transpose(phi_matrix(n,xs))),matrix_transpose(ts)),b);

    bool all_cov = true;
    for (int i = 0; i < prev_mean.size(); i++){
      int a = prev_mean[i][0]*1000, b = mean[i][0]*1000;
      if ( a != b){
        all_cov = false;
      }
    }

    if(all_cov){
      cov_count++;
    }else{
      cov_count = 0;
    }

    if(cov_count == 3){
      cout << ts[0].size() << endl;
      cout << "stop here" << endl;
    }

    prev_mean = mean;

    // cout << "for dubug" << endl;

  }

  return 0;

}