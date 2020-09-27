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

vector<vector<double>> matrix_inverse(vector<vector<double>> matrix){
    vector<vector<double>> ret_matrix;
    return ret_matrix;
}

namespace plt = matplotlibcpp;
int main()
{
    // vector<vector<double>> data = readCSV("raw.csv");
    vector<vector<double>> data = readCSV("test.csv");
    vector<vector<double>> tran_data = matrix_transpose(data);

    // plot raw data
    plt::plot(tran_data.at(0),tran_data.at(1),"ro");
    plt::show();
    vector<vector<double>> ATA = matrix_multiply(tran_data,data);

    vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
}
