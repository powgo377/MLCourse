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

vector<vector<double>> matrix_add(vector<vector<double>> matrixA,vector<vector<double>> matrixB){
    for(int i = 0; i < matrixA.size(); i++){
        for(int j = 0; j < matrixA[0].size(); j++){
            matrixA[i][j] += matrixB[i][j];   
        }
    }
    return matrixA;
}
vector<vector<double>> matrix_minus(vector<vector<double>> matrixA,vector<vector<double>> matrixB){
    for(int i = 0; i < matrixA.size(); i++){
        for(int j = 0; j < matrixA[0].size(); j++){
            matrixA[i][j] -= matrixB[i][j];   
        }
    }
    return matrixA;
}

vector<vector<double>> matrix_multiply_constant(vector<vector<double>> matrix,double constant){
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix[0].size(); j++){
            matrix[i][j] *= constant;   
        }
    }
    return matrix;
}

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
    vector<vector<double>>* LU = LU_decomposition(matrix);
    vector<vector<double>> L = LU [0];
    vector<vector<double>> U = LU [1];

    //construct basis e1, e2, ... and then count yn > xn
    for(int i = 0; i < L.size(); i++){
        vector<double> en(L.size(),0);
        en[i] = 1;
        vector<double> xn = U_linearequation(U,L_linearequation(L,en));
        ret_matrix.push_back(xn);
    }
    
    return matrix_transpose(ret_matrix);
}

vector<double> LSE(int n,double lamda,vector<vector<double>> data){
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


vector<vector<double>> LSE_grediant(vector<vector<double>> A, vector<vector<double>> data,vector<double> startpoint){
    //2ATAx0 - 2ATb
    vector<vector<double>> x0,b;
    x0.push_back(startpoint);
    b.push_back(matrix_transpose(data)[1]);

    vector<vector<double>> ATAx0 = matrix_multiply(matrix_multiply(matrix_transpose(A),A),matrix_transpose(x0));
    vector<vector<double>> ATb = matrix_multiply(matrix_transpose(A),matrix_transpose(b));

    vector<vector<double>> ret = matrix_minus(matrix_multiply_constant(ATAx0,2),matrix_multiply_constant(ATb,2));
    return ret;
}
vector<vector<double>> LSE_Hessian(vector<vector<double>> A){
    //2ATA
    vector<vector<double>> ATA = matrix_multiply(matrix_transpose(A),A);

    // vector<vector<double>> diagATA;
    // for(int i = 0 ; i < ATA.size(); i++){
    //     vector<double> diagATAtmp(A.size(),0);
    //     diagATAtmp[i] = ATA[i][i];
    //     diagATA.push_back(diagATAtmp);
    // }
    // vector<vector<double>> ret = matrix_minus(matrix_multiply_constant(ATA,2),diagATA);
    vector<vector<double>> ret = matrix_multiply_constant(ATA,2);
    return ret;
}
vector<double> Newton(int n,vector<vector<double>> data){
    //start from{1,1,...,1}
    vector<double> argstmp(n,1);
    vector<vector<double>> args;
    args.push_back(argstmp);
    args = matrix_transpose(args);
    vector<vector<double>> A;
    //count matix A
    for(int i = 0; i < data.size(); i++){
        vector<double> row;
        for(int j = 0; j < n; j++){
            row.push_back(pow(data.at(i).at(0) , (n-1-j)));
        }
        A.push_back(row);
    }

    for(int i = 0; i < 10; i++){
        args = matrix_minus(args,matrix_multiply(matrix_inverse_byLU(LSE_Hessian(A)),LSE_grediant(A,data,matrix_transpose(args)[0])));
    }

    return matrix_transpose(args)[0];
}

namespace plt = matplotlibcpp;
int main()
{
    double n = 3,lamda = 10000;
    vector<vector<double>> data = readCSV("raw.csv");
    // vector<vector<double>> data = readCSV("test.csv");
    vector<vector<double>> tran_data = matrix_transpose(data);
    // plot raw data
    plt::plot(tran_data.at(0),tran_data.at(1),"ro");

    //LSE get args
    vector<double> arg_LSE = LSE(n,lamda,data);
    //LSE count loss
    double loss = 0;
    for(int i = 0; i < tran_data[0].size(); i++){
        double perdict = 0;
        for(int j = 0; j < arg_LSE.size(); j++){
            perdict += arg_LSE[j]*pow(tran_data[0][i],arg_LSE.size()-1-j);
        }
        loss += pow(perdict-tran_data[1][i],2);
    }
    //LSE print result
    cout << "LSE:" << endl << "Fitting line:" ;
    for(int j = 0; j < arg_LSE.size(); j++){
            if(j != arg_LSE.size() - 1){
                cout << " " << arg_LSE[j] << "X^" << arg_LSE.size()-1-j << " +";
            } else {
                cout << " " << arg_LSE[j] << endl;
            }
    }
    cout << "Total error: " << loss << endl;
    //LSE plot curve
    double max_x = *max_element(tran_data[0].begin(), tran_data[0].end());
    double min_x = *min_element(tran_data[0].begin(), tran_data[0].end());
    int point_nums = ceil((max_x - min_x)*1.2)+1;
    std::vector<double> x(point_nums),LSE_perdict(point_nums);
    for(int i=0; i < point_nums; i++) {
        x[i] = i-point_nums/2;
        for(int j = 0; j < arg_LSE.size(); j++){
            LSE_perdict[i] += arg_LSE[j]*pow(i-point_nums/2,arg_LSE.size()-1-j);
        }
    }
    plt::plot(x,LSE_perdict);

    //Newton get args
    vector<double> arg_Newton = Newton(n,data);
    //Newton count loss
    loss = 0;
    for(int i = 0; i < tran_data[0].size(); i++){
        double perdict = 0;
        for(int j = 0; j < arg_Newton.size(); j++){
            perdict += arg_Newton[j]*pow(tran_data[0][i],arg_Newton.size()-1-j);
        }
        loss += pow(perdict-tran_data[1][i],2);
    }
    //Newton print result
    cout << "Newton:" << endl << "Fitting line:" ;
    for(int j = 0; j < arg_Newton.size(); j++){
            if(j != arg_Newton.size() - 1){
                cout << " " << arg_Newton[j] << "X^" << arg_Newton.size()-1-j << " +";
            } else {
                cout << " " << arg_Newton[j] << endl;
            }
    }
    cout << "Total error: " << loss << endl;
    //Newton plot curve
    std::vector<double> Newton_perdict(point_nums);
    for(int i=0; i < point_nums; i++) {
        for(int j = 0; j < arg_LSE.size(); j++){
            Newton_perdict[i] += arg_Newton[j]*pow(i-point_nums/2,arg_Newton.size()-1-j);
        }
    }
    plt::plot(x,Newton_perdict);
    
    plt::show();
}
