#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
vector<double> z_sheet{
0.50000,0.50399,0.50798,0.51197,0.51595,0.51994,0.52392,0.52790,0.53188,0.53586,
0.53983,0.54380,0.54776,0.55172,0.55567,0.55962,0.56356,0.56749,0.57142,0.57535,
0.57926,0.58317,0.58706,0.59095,0.59483,0.59871,0.60257,0.60642,0.61026,0.61409,
0.61791,0.62172,0.62552,0.62930,0.63307,0.63683,0.64058,0.64431,0.64803,0.65173,
0.65542,0.65910,0.66276,0.66640,0.67003,0.67364,0.67724,0.68082,0.68439,0.68793,
0.69146,0.69497,0.69847,0.70194,0.70540,0.70884,0.71226,0.71566,0.71904,0.72240,
0.72575,0.72907,0.73237,0.73565,0.73891,0.74215,0.74537,0.74857,0.75175,0.75490,
0.75804,0.76115,0.76424,0.76730,0.77035,0.77337,0.77637,0.77935,0.78230,0.78524,
0.78814,0.79103,0.79389,0.79673,0.79955,0.80234,0.80511,0.80785,0.81057,0.81327,
0.81594,0.81859,0.82121,0.82381,0.82639,0.82894,0.83147,0.83398,0.83646,0.83891,
0.84134,0.84375,0.84614,0.84849,0.85083,0.85314,0.85543,0.85769,0.85993,0.86214,
0.86433,0.86650,0.86864,0.87076,0.87286,0.87493,0.87698,0.87900,0.88100,0.88298,
0.88493,0.88686,0.88877,0.89065,0.89251,0.89435,0.89617,0.89796,0.89973,0.90147,
0.90320,0.90490,0.90658,0.90824,0.90988,0.91149,0.91309,0.91466,0.91621,0.91774,
0.91924,0.92073,0.92220,0.92364,0.92507,0.92647,0.92785,0.92922,0.93056,0.93189,
0.93319,0.93448,0.93574,0.93699,0.93822,0.93943,0.94062,0.94179,0.94295,0.94408,
0.94520,0.94630,0.94738,0.94845,0.94950,0.95053,0.95154,0.95254,0.95352,0.95449,
0.95543,0.95637,0.95728,0.95818,0.95907,0.95994,0.96080,0.96164,0.96246,0.96327,
0.96407,0.96485,0.96562,0.96638,0.96712,0.96784,0.96856,0.96926,0.96995,0.97062,
0.97128,0.97193,0.97257,0.97320,0.97381,0.97441,0.97500,0.97558,0.97615,0.97670,
0.97725,0.97778,0.97831,0.97882,0.97932,0.97982,0.98030,0.98077,0.98124,0.98169,
0.98214,0.98257,0.98300,0.98341,0.98382,0.98422,0.98461,0.98500,0.98537,0.98574,
0.98610,0.98645,0.98679,0.98713,0.98745,0.98778,0.98809,0.98840,0.98870,0.98899,
0.98928,0.98956,0.98983,0.99010,0.99036,0.99061,0.99086,0.99111,0.99134,0.99158,
0.99180,0.99202,0.99224,0.99245,0.99266,0.99286,0.99305,0.99324,0.99343,0.99361,
0.99379,0.99396,0.99413,0.99430,0.99446,0.99461,0.99477,0.99492,0.99506,0.99520,
0.99534,0.99547,0.99560,0.99573,0.99585,0.99598,0.99609,0.99621,0.99632,0.99643,
0.99653,0.99664,0.99674,0.99683,0.99693,0.99702,0.99711,0.99720,0.99728,0.99736,
0.99744,0.99752,0.99760,0.99767,0.99774,0.99781,0.99788,0.99795,0.99801,0.99807,
0.99813,0.99819,0.99825,0.99831,0.99836,0.99841,0.99846,0.99851,0.99856,0.99861,
0.99865,0.99869,0.99874,0.99878,0.99882,0.99886,0.99889,0.99893,0.99896,0.99900,
0.99903,0.99906,0.99910,0.99913,0.99916,0.99918,0.99921,0.99924,0.99926,0.99929,
0.99931,0.99934,0.99936,0.99938,0.99940,0.99942,0.99944,0.99946,0.99948,0.99950,
0.99952,0.99953,0.99955,0.99957,0.99958,0.99960,0.99961,0.99962,0.99964,0.99965,
0.99966,0.99968,0.99969,0.99970,0.99971,0.99972,0.99973,0.99974,0.99975,0.99976,
0.99977,0.99978,0.99978,0.99979,0.99980,0.99981,0.99981,0.99982,0.99983,0.99983,
0.99984,0.99985,0.99985,0.99986,0.99986,0.99987,0.99987,0.99988,0.99988,0.99989,
0.99989,0.99990,0.99990,0.99990,0.99991,0.99991,0.99992,0.99992,0.99992,0.99992,
0.99993,0.99993,0.99993,0.99994,0.99994,0.99994,0.99994,0.99995,0.99995,0.99995,
0.99995,0.99995,0.99996,0.99996,0.99996,0.99996,0.99996,0.99996,0.99997,0.99997};

int pic_num = 0;
vector<int> num_count(10,0);
int row_num = 0;
int column_num = 0;
int read_Integer(ifstream & in){
    int result = 0;
    char readbyte;
    for(int i = 0; i < 4; i++){
        result *= 256;
        in.get(readbyte);
        int tmp = int(readbyte);
        if(readbyte < 0){    
            result += tmp + 256;
        }else{
            result += tmp;
        }
    }
    return result;
}
vector<vector<vector<char>>> read_binary_to_pics(string image_file){
    vector<vector<vector<char>>> ret;
    char readbyte;
    ifstream in(image_file , ios::in | ios::binary);

    int magic_num = read_Integer(in);
    pic_num = read_Integer(in);
    row_num = read_Integer(in);
    column_num = read_Integer(in);

    for(int i = 0; i < pic_num ;i++){
        vector<vector<char>> tmp_picture;
        for(int j = 0; j < row_num ; j++){
            vector<char> tmp_row;
            for(int k = 0; k < column_num ; k++){
                in.get(readbyte);
                tmp_row.push_back(readbyte);
            }
            tmp_picture.push_back(tmp_row);
        }  
        ret.push_back(tmp_picture);
    }
    int size = ret.size();
    return ret;
}
vector<char> read_binary_to_labels(string image_file){
    vector<char> ret;
    char readbyte;
    ifstream in(image_file , ios::in | ios::binary);

    int magic_num = read_Integer(in);
    int label_num = read_Integer(in);
    
    for(int i = 0; i < label_num ;i++){
        in.get(readbyte);
        num_count[int(readbyte)] += 1;
        ret.push_back(readbyte); 
    }
    int size = ret.size();
    return ret;
}
vector<vector<vector<char>>> read_binary_to_t_pics(string image_file){
    vector<vector<vector<char>>> ret;
    char readbyte;
    ifstream in(image_file , ios::in | ios::binary);

    int magic_num = read_Integer(in);
    int pic_num = read_Integer(in);
    int row_num = read_Integer(in);
    int column_num = read_Integer(in);

    for(int i = 0; i < pic_num ;i++){
        vector<vector<char>> tmp_picture;
        for(int j = 0; j < row_num ; j++){
            vector<char> tmp_row;
            for(int k = 0; k < column_num ; k++){
                in.get(readbyte);
                tmp_row.push_back(readbyte);
            }
            tmp_picture.push_back(tmp_row);
        }  
        ret.push_back(tmp_picture);
    }
    int size = ret.size();
    return ret;
}
vector<char> read_binary_to_t_labels(string image_file){
    vector<char> ret;
    char readbyte;
    ifstream in(image_file , ios::in | ios::binary);

    int magic_num = read_Integer(in);
    int label_num = read_Integer(in);
    
    for(int i = 0; i < label_num ;i++){
        in.get(readbyte);
        ret.push_back(readbyte); 
    }
    int size = ret.size();
    return ret;
}

//int/row/col/bin
vector<vector<vector<vector<double>>>> pre_process_data(vector<vector<vector<char>>> &pictures, vector<char> &labels,bool discrete){
    //// init 0s to /int/row/col/bin
    vector<double> bins(32, 0);
    vector<vector<double>> cols(column_num,bins);
    vector<vector<vector<double>>> rows(row_num,cols);
    vector<vector<vector<vector<double>>>> ret(10,rows);

    //// count into bins
    for(int i = 0; i < labels.size(); i++){
        int label = labels[i];
        for(int j = 0; j < row_num; j++){
            // cout << "1234" << endl;
            for(int k = 0; k < column_num; k++){
                // cout << "1234" << endl;
                int bin_rank = 0;
                if(int(pictures[i][j][k])< 0){
                    bin_rank = (int(pictures[i][j][k])+256)/8;
                }else{
                    bin_rank = (int(pictures[i][j][k]))/8;    
                }
                ret[label][j][k][bin_rank] += 1;
            }
        }
    }

    pictures.clear();
    vector<vector<vector<char>>>().swap(pictures);
    labels.clear();
    vector<char>().swap(labels);

    //// set 0 bin to min value
    // find min
    if(discrete){
        int min = pic_num;
        for(int label = 0; label < 10; label++){
            for(int i = 0; i < row_num; i++){
                for(int j = 0; j < column_num; j++){
                    for(int k = 0; k < 32; k++){
                        int tmp = ret[label][i][j][k];
                        if( tmp > 0 && tmp < min){
                            min = tmp;
                        }
                    }
                }
            }
        }
        
        // set min
        for(int label = 0; label < 10; label++){
            for(int i = 0; i < row_num; i++){
                for(int j = 0; j < column_num; j++){
                    for(int k = 0; k < 32; k++){
                        if( ret[label][i][j][k] == 0 ){
                            ret[label][i][j][k] = min;
                        }
                    }
                }
            }
        }
    }    
    return ret;
}
void pre_process_data_gaussian(vector<vector<vector<vector<double>>>> &pre){
    for(int label =0; label < 10; label++){
        for(int i = 0; i < row_num; i++){
            for(int j = 0; j < column_num; j++){
                int count = 0;
                double mean = 0;
                long double deviation = 0;
                for(int k = 0; k < 32; k++){
                    count += pre[label][i][j][k];
                    mean += pre[label][i][j][k] * k ;
                }
                mean = mean / count;
                for(int k = 0; k < 32; k++){
                    deviation += pre[label][i][j][k]*(k-mean)*(k-mean);
                }
                deviation = deviation / count;
                deviation = sqrt(deviation);

                double sum = 0;

                for(int k = 0; k < 32; k++){
                    double prob;
                    ////查Z表法
                    // if(k<mean){
                    // double x = mean+ mean - k + 0.5;
                    if(pre[label][i][j][k] != 0 && deviation == 0){
                        prob = 1;
                    }else if (pre[label][i][j][k] == 0 && deviation == 0)
                    {
                        prob = 0;
                    }else{
                        double probm5 = 0,probp5 = 0;
                        double zm05 = (k-0.5 - mean)/deviation;
                        if (zm05 >= -4.0 && zm05 <= 4.0){
                            int z_index = abs((zm05*100)/1);
                            if(k-0.5<mean){
                                probm5 = 1 - z_sheet[z_index];
                            }else{
                                probm5 = z_sheet[z_index];
                            }
                        }else if(zm05>4.0){
                            probm5 = 1;
                        }
                        double zp05 = (k+0.5 - mean)/deviation;
                        if (zp05 >= -4.0 && zp05 <= 4.0){
                            int z_index = abs((zp05*100)/1);
                            if(k+0.5<mean){
                                probp5 = 1 - z_sheet[z_index];
                            }else{
                                probp5 = z_sheet[z_index];
                            }
                        }else if(zp05>4.0){
                            probp5 = 1;
                        }
                        prob = probp5-probm5;
                        pre[label][i][j][k] = prob*count;
                    }
                }
                // cout << "456" << endl;
            }
        }
    }
    double min = pic_num;
    for(int label = 0; label < 10; label++){
        for(int i = 0; i < row_num; i++){
            for(int j = 0; j < column_num; j++){
                for(int k = 0; k < 32; k++){
                    double tmp = pre[label][i][j][k];
                    if( tmp > 0 && tmp < min){
                        min = tmp;
                    }
                }
            }
        }
    }
    
    // set min
    for(int label = 0; label < 10; label++){
        for(int i = 0; i < row_num; i++){
            for(int j = 0; j < column_num; j++){
                for(int k = 0; k < 32; k++){
                    if( pre[label][i][j][k] == 0 ){
                        pre[label][i][j][k] = min;
                    }
                }
            }
        }
    }
}

double classfied(vector<vector<vector<vector<double>>>> &pic_classified_in_bins,vector<vector<vector<char>>> &t_pictures,vector<char> &t_labels){
    int error_count = 0;
    for(int i = 0; i < t_labels.size(); i++){
        vector<long double> likelis(10,0);
        for(int j = 0; j < 10; j++){
            long double likeli = num_count[j]/static_cast<long double>(pic_num);
            for(int k = 0; k < 28; k++){
                // cout << "tmp" << endl;
                for(int l = 0; l < 28; l++){
                    // cout << "tmp" << endl;
                    int bin_level = 0;
                    if(int(t_pictures[i][k][l])< 0){
                        bin_level = (int(t_pictures[i][k][l])+256)/8;
                    }else{
                        bin_level = (int(t_pictures[i][k][l]))/8;    
                    }
                    double bin_p = pic_classified_in_bins[j][k][l][bin_level];
                    long double j_count = static_cast<long double>(num_count[j]);
                    likeli *= bin_p/j_count;
                }
            }
            likelis[j] = likeli;
        }
        cout << "Postirior (in log scale):" << endl;
        vector<double> log_scale(0,0);
        long double sum = 0;
        for(int i = 0; i < 10; i++){
            log_scale.push_back(log(likelis[i])*(-1));
            long double tmp1 = likelis[i];
            long double tmp2 = log(likelis[i]);
            sum += log(likelis[i])*(-1);
        }
        double min = sum;
        int min_count = 10;
        for(int i = 0; i < 10; i++){
            log_scale[i] = log_scale[i] / sum;
            cout << i << ": " << log_scale[i] << endl;
            if(log_scale[i] < min){
                min_count = i;
                min = log_scale[i];
            }
        }
        cout << "Prediction: "<< min_count <<", Ans: "<<int(t_labels[i])<< endl << endl;
        if(min_count != int(t_labels[i])){
            error_count += 1;
        }
    }
    return error_count/static_cast<long double>(t_labels.size());
}
vector<vector<vector<vector<double>>>> print_num(vector<vector<vector<vector<double>>>> &bins, bool discrete){
    //// count into bins
    for(int i = 0; i < 10; i++){
        cout << "   " << endl;
        for(int j = 0; j < row_num; j++){
            for(int k = 0; k < column_num; k++){
                if(discrete){
                    int count = 0;
                    for(int l = 0; l < 16; l++){
                        count += bins[i][j][k][32-l-1];
                    }

                    if(count >= num_count[i]/2){
                        cout << "1";
                    }else{
                        cout << "0";
                    }
                }else{
                    int max = 0,max_index = 0;
                    bool is_black = false;
                    for(int l = 0; l < 32; l++){
                        if (bins[i][j][k][l] > max){
                            max = bins[i][j][k][l];
                            max_index = l;
                        }
                    }
                    if(max_index>15){
                        is_black = true;
                    }

                    if(is_black){
                        cout << "1";
                    }else{
                        cout << "0";
                    }
                }
            }
            cout <<" :"<< j << endl;
        }
    }
}
int main()
{
    bool is_discrete = false;
    vector<vector<vector<char>>> pictures = read_binary_to_pics("train-images.idx3-ubyte");
    vector<char> labels = read_binary_to_labels("train-labels.idx1-ubyte");
    
    vector<vector<vector<vector<double>>>> pic_classified_in_bins = pre_process_data(pictures,labels,is_discrete);

    ////change to continuous
    if(!is_discrete){
        pre_process_data_gaussian(pic_classified_in_bins);
    }
    
    vector<vector<vector<char>>> t_pictures = read_binary_to_t_pics("t10k-images.idx3-ubyte");
    vector<char> t_labels = read_binary_to_t_labels("t10k-labels.idx1-ubyte");
    
    double error_rate = classfied(pic_classified_in_bins,t_pictures,t_labels);
    print_num(pic_classified_in_bins,false);

    cout << "Error rate: "<< error_rate << endl;

    return 0;
}
