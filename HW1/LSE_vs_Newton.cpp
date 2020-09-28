#include "../../lib/matplotlib-cpp-master/matplotlibcpp.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

vector<vector<double>> readCSV(string csvfile){
    vector<vector<double>> data;
    
    string line;
    fstream fileStream;
    fileStream.open(csvfile);
    
    while (getline( fileStream, line,'\n')) 
	{
        vector<double> datapoint;
        string value;
        stringstream lineStream(line); 

        while(getline(lineStream, value, ',')){
            //TODO:stod get different value
            datapoint.push_back(stod(value));
            // cout << setprecision(100)<<datapoint.back()<<',';
        }
        data.push_back(datapoint);
        // cout << endl;
	}
    return data;
}

vector<vector<double>> matrix_transpose(vector<vector<double>> matrix){
    vector<vector<double>> tranposed_matrix;
    for(int i = 0; i < matrix.at(0).size(); i++){
        vector<double> row;
        for(int j = 0; j < matrix.size(); j++){
            row.push_back(matrix.at(j).at(i));
        }
        tranposed_matrix.push_back(row);
    }
    return tranposed_matrix;
}
////square mult square
vector<vector<double>> matrix_multiply(vector<vector<double>> matrixA,vector<vector<double>> matrixB){
    vector<vector<double>> ret_matrix;
    for(int i = 0; i < matrixA.size(); i++){
        vector<double> row;
        for(int j = 0; j < matrixB.at(0).size(); j++){
            double mult = 0;
            for(int k = 0; k < matrixB.size(); k++){
                mult += matrixA.at(i).at(k) * matrixB.at(k).at(j);
            }
            row.push_back(mult);
        }
        ret_matrix.push_back(row);
    }
    return ret_matrix;
}

vector<vector<double>> LU[2];
vector<vector<double>>* LU_decomposition(vector<vector<double>> matrix){
    vector<vector<double>> L;
    vector<vector<double>> U = matrix;
    for(int i = 0; i < matrix.size(); i++){
        vector<double> L_row(matrix.size(),0);
        L_row[i] = 1;
        for(int j = 0; j < i; j++){
            // cout << "LU counter : " << i << j << endl;
            double ratio = U[i][j]/U[j][j];
            U[i][j] = 0;
            for(int k = j+1; k < matrix.size(); k++){
                U[i][k] = U[i][k] - ratio*U[j][k];
            }
            L_row[j] = ratio;
        }
        L.push_back(L_row);
    }
    LU[0] = L;
    LU[1] = U;
    return LU;
}

vector<double> L_linearequation(vector<vector<double>> matrix, vector<double> b){
    for(int i = 0; i < b.size(); i++){       
        for(int j = 0; j < i; j++){
            b[i] = b[i] - matrix[i][j]*b[j];
            matrix[i][j] = 0;
        }
        b[i] = b[i] / matrix[i][i];
        matrix[i][i] = 1;
    }
    return b;
}
vector<double> U_linearequation(vector<vector<double>> matrix, vector<double> b){
    for(int i = 0; i < b.size(); i++){       
        for(int j = 0; j < i; j++){
            b[b.size()-1-i] = b[b.size()-1-i] - matrix[b.size()-1-i][b.size()-1-j]*b[b.size()-1-j];
            matrix[b.size()-1-i][b.size()-1-j] = 0;
        }
        b[b.size()-1-i] = b[b.size()-1-i] / matrix[b.size()-1-i][b.size()-1-i];
        matrix[b.size()-1-i][b.size()-1-i] = 1;
    }
    return b;
}

vector<vector<double>> matrix_inverse_byLU(vector<vector<double>> matrix){
    vector<vector<double>> ret_matrix;
    // vector<double> row1 = {3,-1,2};
    // vector<double> row2 = {6,-1,5};
    // vector<double> row3 = {-9,7,3};
    // // vector<double> b = {10,2,15};
    // vector<vector<double>> matrix1;
    // matrix1.push_back(row1);
    // matrix1.push_back(row2);
    // matrix1.push_back(row3);
    vector<vector<double>>* LU = LU_decomposition(matrix);
    vector<vector<double>> L = LU [0];
    vector<vector<double>> U = LU [1];

    vector<double> e1 = {1,0,0};
    vector<double> e2 = {0,1,0};
    vector<double> e3 = {0,0,1};

    vector<double> y1 = L_linearequation(L,e1);
    vector<double> y2 = L_linearequation(L,e2);
    vector<double> y3 = L_linearequation(L,e3);
    vector<double> x1 = U_linearequation(U,y1);
    vector<double> x2 = U_linearequation(U,y2);
    vector<double> x3 = U_linearequation(U,y3);
    
    ret_matrix.push_back(x1);
    ret_matrix.push_back(x2);
    ret_matrix.push_back(x3);

    vector<vector<double>> tmp = matrix_multiply(matrix,matrix_transpose(ret_matrix));

    return matrix_transpose(ret_matrix);
}

vector<double> LSE(int n,double lamda,vector<vector<double>> data){
    vector<double> args;

    vector<vector<double>> A;
    //count matix A
    for(int i = 0; i < data.size(); i++){
        vector<double> row;
        for(int j = 0; j < n; j++){
            row.push_back(pow(data.at(i).at(0) , (n-1-j)));
        }
        A.push_back(row);
    }

    vector<vector<double>> tran_A = matrix_transpose(A);
    vector<vector<double>> ATA = matrix_multiply(tran_A,A);

    //ATA plus lamda*I
    for(int i = 0; i < ATA.size(); i++){
        ATA[i][i] = ATA[i][i] + lamda;
    }

    vector<vector<double>> ATA_plusL_inverse = matrix_inverse_byLU(ATA);
    vector<vector<double>> ATA_plusL_inverseAT,b;
    b.push_back(matrix_transpose(data)[1]);
    ATA_plusL_inverseAT = matrix_multiply(ATA_plusL_inverse,matrix_transpose(A));
    vector<vector<double>> tmp = matrix_multiply(ATA_plusL_inverseAT,matrix_transpose(b));

    return matrix_transpose(tmp)[0];
}

namespace plt = matplotlibcpp;
int main()
{
    vector<vector<double>> data = readCSV("raw.csv");
    // vector<vector<double>> data = readCSV("test.csv");

    vector<double> arg_LSE = LSE(3,10000,data);

    // plot raw data
    vector<vector<double>> tran_data = matrix_transpose(data);
    plt::plot(tran_data.at(0),tran_data.at(1),"ro");
    // plt::show();
}
