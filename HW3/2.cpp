#include <iostream>
#include <random>

using namespace std;

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

int main()
{
  //////// 1a
  ////input 1a
  double mean = 50,deviation = 10;

  cout << "Mean : " << mean;
  cout << "  Deviation : " << deviation;
  cout << "  Desult : " << gaussian_rand(mean,deviation) << endl;

  //////// 1b
  ////input 1b
  int n = 2;
  double a = 1;
  vector<double> w{2,3};

  mt19937 rng(std::random_device{}());
  uniform_real_distribution<double> x_distribution(-1.0,1.0);
  double x = x_distribution(rng);
  double y = 0;
  for(int i = 0; i < n; i++){
    y += w[i]*pow(x,i);
  }
  y += gaussian_rand(0,a);
  cout << "Y : " << y << endl;

  //////// 2
  double x1 = gaussian_rand(mean,deviation);
  int x_count = 1,cov_count = 0;
  double mean_guess = x1, deviation_guess = 0,s = 0;
  int mean_cov = 0,deviation_cov = 0;

  while(true && cov_count < 3){
    double xn = gaussian_rand(mean,deviation);
    double meann = ( mean_guess * x_count + xn )/( x_count + 1 );
    s = s + ( xn - mean_guess )*( xn - meann );

    mean_guess = meann;
    deviation_guess = sqrt(s/x_count+1);
    x_count++; 

    int mean_cov_tmp = mean_guess*10000;
    int deviation_cov_tmp = deviation_guess*10000;
    if(mean_cov_tmp == mean_cov && deviation_cov_tmp == deviation_cov){
      cov_count++;
    }else{
      cov_count = 0;
      mean_cov = mean_cov_tmp;
      deviation_cov = deviation_cov_tmp;
    }
  }
  cout << "mean_guess : " << mean_guess <<" deviation_guess : "<< deviation_guess <<" iteration count : " << x_count << endl;

  return 0;
}