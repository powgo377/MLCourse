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

namespace plt = matplotlibcpp;
int main()
{
    vector<vector<double>> data = readCSV("raw.csv");
    vector<double> x_array;
    vector<double> y_array;

    for (int i = 0; i < data.size(); i++){
        x_array.push_back(data.at(i).at(0));
        y_array.push_back(data.at(i).at(1));
        cout << data.at(i).at(0) << ',' << data.at(i).at(1) << endl;
    }
    plt::plot(x_array,y_array,"ro");
    plt::show();

    vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
}
