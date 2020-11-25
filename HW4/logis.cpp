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

vector<vector<double>> sigmoid(vector<vector<double>> X, vector<vector<double>> W){
    vector<vector<double>> ret;
    vector<vector<double>> XW = matrix_multiply(X, W);
    for (int i = 0; i < XW.size(); i++){
        vector<double> tmp;
        tmp.push_back(1.0 / (1.0 + exp(-1.0 * XW[i][0])));
        ret.push_back(tmp);
    }
    return ret;
}

namespace plt = matplotlibcpp;
int main()
{
    // double N = 50;
    // double mx1 = 1;
    // double my1 = 1;
    // double mx2 = 10;
    // double my2 = 10;
    // double vx1 = 2;
    // double vy1 = 2;
    // double vx2 = 2;
    // double vy2 = 2;
    
    double N = 50;
    double mx1 = 1;
    double my1 = 1;
    double mx2 = 3;
    double my2 = 3;
    double vx1 = 2;
    double vy1 = 2;
    double vx2 = 4;
    double vy2 = 4;

	vector<double> empty;
	vector<vector<double>> x1s;
	vector<vector<double>> y1s;
	vector<vector<double>> x2s;
	vector<vector<double>> y2s;
    vector<vector<double>> X;
    vector<vector<double>> Y;
	x1s.push_back(empty);
	y1s.push_back(empty);
	x2s.push_back(empty);
	y2s.push_back(empty);
    // cout << "Mean : " << mx1;
    // cout << "  Deviation : " << sqrt(vx1);
    for(int i = 0;i < N;i++){
        double x1 = gaussian_rand(mx1,vx1);
        x1s[0].push_back(x1);
        double y1 = gaussian_rand(my1,vy1);
        y1s[0].push_back(y1);
        vector<double> tmp;
        tmp.push_back(x1);
        tmp.push_back(y1);
        tmp.push_back(1);
        X.push_back(tmp);
        vector<double> tmpy0(1,0);
        Y.push_back(tmpy0);
        double x2 = gaussian_rand(mx2,vx2);
        x2s[0].push_back(x2);
        double y2 = gaussian_rand(my2,vy2);
        y2s[0].push_back(y2);
        vector<double> tmp2;
        tmp2.push_back(x2);
        tmp2.push_back(y2);
        tmp2.push_back(1);
        X.push_back(tmp2);
        vector<double> tmpy1(1,1);
        Y.push_back(tmpy1);
    }

    vector<double> W_row(1,1);
    vector<vector<double>> W(3,W_row);
    vector<vector<double>> prev_W = W;
    
    ////gradient version
    // int count = 0;
    int cov_count = 0;
    bool not_diff=true;
    while(not_diff){
        // count++;s
        vector<vector<double>> sigmoid_XW = sigmoid(X,W);
        vector<vector<double>> gradient_w = matrix_multiply(matrix_transpose(X),matrix_minus(Y,sigmoid_XW));
        W = matrix_add(W,matrix_multiply_constant(gradient_w,0.001));
        if((int(W[0][0]*1000) == int(prev_W[0][0]*1000)) && ((int(W[1][0]*1000) == int(prev_W[1][0]*1000))) && ((int(W[2][0]*1000) == int(prev_W[2][0]*1000)))){
            cov_count++;
        }else{
            cov_count = 0;
        }
        if(cov_count == 3){
            not_diff = false;
        }
        prev_W = W;
    }


    vector<vector<double>> W_new(3,W_row);
    vector<vector<double>> prev_W_new = W_new;
    vector<vector<double>> pprev_W_new = W_new;
    ////newton's version
    int count = 0;
    int cov_count_new = 0;
    bool not_diff_new=true;
    bool boom = false;
    while(not_diff_new){
        count++;
        vector<vector<double>> D;
        for(int i = 0; i < X.size(); i++){
            vector<double> tmp(X.size(),0);
            double a =-1*(X[i][0]*W_new[0][0] + X[i][1]*W_new[1][0] + X[i][2]*W_new[2][0]); 
            if(a>500){
                a = 500;
            }
            double exp_mXW = exp(a);
            tmp[i] = exp_mXW / pow((1 + exp_mXW),2);
            D.push_back(tmp);
        }
        vector<vector<double>> H = matrix_multiply(matrix_transpose(X), matrix_multiply(D, X));
        vector<vector<double>> gradient_w = matrix_multiply(matrix_transpose(X),matrix_minus(Y,sigmoid(X,W_new)));
        vector<vector<double>> H_inv = matrix_inverse_byLU(H);
        // double det = H[0][0]*H[1][1]-H[1][0]*H[0][1];
        double det = H[0][0] * (H[1][1] * H[2][2] - H[2][1] * H[1][2]) - H[0][1] * (H[1][0] * H[2][2] - H[1][2] * H[2][0]) + H[0][2] * (H[1][0] * H[2][1] - H[1][1] * H[2][0]);
        int test = 0;
        if(boom || abs(det)<0.00000001){
            ////steepest gdescent
            W_new = matrix_add(W_new,matrix_multiply_constant(gradient_w,0.001));
        }else{
            ////Newtons
            vector<vector<double>> move = matrix_multiply(matrix_inverse_byLU(H), gradient_w);
            if(abs(move[0][0])>200){
                // W_new = matrix_add(W_new,matrix_multiply_constant(gradient_w,0.01));
                boom = true;
            }else{
                W_new = matrix_add(W_new,move);
            }
        }
        
        if((int(W_new[0][0]*10) == int(prev_W_new[0][0]*10)) && ((int(W_new[1][0]*10) == int(prev_W_new[1][0]*10))) && ((int(W_new[2][0]*10) == int(prev_W_new[2][0]*10)))){
            cov_count_new++;
        }else{
            cov_count_new = 0;
        }
        
        if(cov_count_new == 3){
            not_diff_new = false;
        }
        pprev_W_new = prev_W_new;
        prev_W_new = W_new;
    }
    


    ////plot
    plt::subplot(1,3,1);
    plt::plot(x1s[0],y1s[0],"ro");
    plt::plot(x2s[0],y2s[0],"bo");
    
    ////result_gradient
    double p11=0,p12=0,p21=0,p22=0;
    cout << "Gradient descent:" << endl << endl << "w:" << endl;
    cout << W[0][0] << endl << W[1][0] << endl << W[2][0] << endl << endl;
    plt::subplot(1,3,2);
    vector<double> rx1,ry1,rx2,ry2;
    vector<vector<double>> result_Gradient = sigmoid(X, W);
    for(int i = 0;i < X.size(); i++){
        if(result_Gradient[i][0] < 0.5){
            if(Y[i][0] == 0){
                p11++;
            }else{
                p12++;
            }
            rx1.push_back(X[i][0]);
            ry1.push_back(X[i][1]);
        }else{
            if(Y[i][0] == 0){
                p21++;
            }else{
                p22++;
            }
            rx2.push_back(X[i][0]);
            ry2.push_back(X[i][1]);
        }
    }
    cout << "Confusion Matrix:" << endl << "            Predict cluster 1 Predict cluster 2" << endl;
    cout << "Is cluster 1       "<< p11 << "             " << p21 << endl;
    cout << "Is cluster 2       "<< p12 << "             " << p22 << endl << endl;
    cout << "Sensitivity (Successfully predict cluster 1):" << p11/(p11+p12) << endl;
    cout << "Sensitivity (Successfully predict cluster 2):" << p22/(p21+p22) << endl << endl;
    cout << "---------------------------------------------------------" << endl;
    plt::plot(rx1,ry1,"ro");
    plt::plot(rx2,ry2,"bo");

    ////plot_newtons
    p11=0,p12=0,p21=0,p22=0;
    cout << "Newton's method:" << endl << endl << "w:" << endl;
    cout << W_new[0][0] << endl << W_new[1][0] << endl << W_new[2][0] << endl << endl;
    plt::subplot(1,3,3);
    vector<double> nrx1,nry1,nrx2,nry2;
    vector<vector<double>> nresult_Gradient = sigmoid(X, W_new);
    for(int i = 0;i < X.size(); i++){
        if(nresult_Gradient[i][0] < 0.5){
            if(Y[i][0] == 0){
                p11++;
            }else{
                p12++;
            }
            nrx1.push_back(X[i][0]);
            nry1.push_back(X[i][1]);
        }else{
            if(Y[i][0] == 0){
                p21++;
            }else{
                p22++;
            }
            nrx2.push_back(X[i][0]);
            nry2.push_back(X[i][1]);
        }
    }
    cout << "Confusion Matrix:" << endl << "            Predict cluster 1 Predict cluster 2" << endl;
    cout << "Is cluster 1       "<< p11 << "             " << p21 << endl;
    cout << "Is cluster 2       "<< p12 << "             " << p22 << endl << endl;
    cout << "Sensitivity (Successfully predict cluster 1):" << p11/(p11+p12) << endl;
    cout << "Sensitivity (Successfully predict cluster 2):" << p22/(p21+p22) << endl;
    plt::plot(nrx1,nry1,"ro");
    plt::plot(nrx2,nry2,"bo");
    plt::show();

    return 0;
}