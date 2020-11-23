#include <iostream>
#include <random>
#include "./LSE_vs_Newton.cpp"
#include "../../lib/matplotlib-cpp-master/matplotlibcpp.h"

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

namespace plt = matplotlibcpp;
int main()
{
    double N = 50;
    double mx1 = 1;
    double mx2 = 10;
    double my1 = 1;
    double my2 = 10;
    double vx1 = 2;
    double vx2 = 2;
    double vy1 = 2;
    double vy2 = 2;
    
	vector<double> empty;
	vector<vector<double>> x1s;
	vector<vector<double>> y1s;
	vector<vector<double>> x2s;
	vector<vector<double>> y2s;
	x1s.push_back(empty);
	y1s.push_back(empty);
	x2s.push_back(empty);
	y2s.push_back(empty);
    // cout << "Mean : " << mx1;
    // cout << "  Deviation : " << sqrt(vx1);
    for(int i = 0;i < N;i++){
		x1s[0].push_back(gaussian_rand(mx1,vx1));
		y1s[0].push_back(gaussian_rand(my1,vy1));
		x2s[0].push_back(gaussian_rand(mx2,vx2));
		y2s[0].push_back(gaussian_rand(my2,vy2));
    }
	plt::plot(x1s[0],y1s[0],"ro");
	plt::plot(x2s[0],y2s[0],"bo");
	plt::show();
    cout << "  Result : " << gaussian_rand(mx1,vx1) << endl;

    return 0;
}