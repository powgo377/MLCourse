#include "../../lib/matplotlib-cpp-master/matplotlibcpp.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

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
    // int count = 0;
    // while(!in.eof()) {
    //     in.get(readbyte);
    //     count += 1;
    // }
    // cout << count << ":"<<int(readbyte) << endl;
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
    // int count = 0;
    // while(!in.eof()) {
    //     in.get(readbyte);
    //     count += 1;
    // }
    // cout << count << ":"<<int(readbyte) << endl;
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
    // int count = 0;
    // while(!in.eof()) {
    //     in.get(readbyte);
    //     count += 1;
    // }
    // cout << count << ":"<<int(readbyte) << endl;
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
    // int count = 0;
    // while(!in.eof()) {
    //     in.get(readbyte);
    //     count += 1;
    // }
    // cout << count << ":"<<int(readbyte) << endl;
    return ret;
}

//int/row/col/bin
vector<vector<vector<vector<int>>>> pre_process_data(vector<vector<vector<char>>> &pictures, vector<char> &labels){
    //// init 0s to /int/row/col/bin
    vector<int> bins(32, 0);
    vector<vector<int>> cols(column_num,bins);
    vector<vector<vector<int>>> rows(row_num,cols);
    vector<vector<vector<vector<int>>>> ret(10,rows);

    //// count into bins
    for(int i = 0; i < labels.size(); i++){
        int label = labels[i];
        for(int j = 0; j < row_num; j++){
            // cout << "1234" << endl;
            for(int k = 0; k < column_num; k++){
                // cout << "1234" << endl;
                int bin_rank = 0;
                if(int(pictures[i][j][k])< 0){
                    bin_rank = (int(pictures[i][j][k])+256)/32;
                }else{
                    bin_rank = (int(pictures[i][j][k]))/32;    
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
    return ret;
}

void discrete_test(vector<vector<vector<vector<int>>>> &pic_classified_in_bins,vector<vector<vector<char>>> &t_pictures,vector<char> &t_labels){
    for(int i = 0; i < t_labels.size(); i++){
        vector<long double> likelis(10,0);
        for(int j = 0; j < 10; j++){
            long double likeli = num_count[j]/static_cast<long double>(pic_num);
            for(int k = 0; k < 28; k++){
                for(int l = 0; l < 28; l++){
                    int bin_level = 0;
                    if(int(t_pictures[i][k][l])< 0){
                        bin_level = (int(t_pictures[i][k][l])+256)/32;
                    }else{
                        bin_level = (int(t_pictures[i][k][l]))/32;    
                    }
                    int bin_p = pic_classified_in_bins[j][k][l][bin_level];
                    long double j_count = static_cast<long double>(num_count[j]);
                    likeli *= bin_p/j_count;
                }
            }
            likelis[j] = likeli;
        }
        cout << "1234" << endl;
    }
}
int main()
{
    vector<vector<vector<char>>> pictures = read_binary_to_pics("train-images.idx3-ubyte");
    vector<char> labels = read_binary_to_labels("train-labels.idx1-ubyte");
    
    vector<vector<vector<vector<int>>>> pic_classified_in_bins = pre_process_data(pictures,labels);

    vector<vector<vector<char>>> t_pictures = read_binary_to_t_pics("t10k-images.idx3-ubyte");
    vector<char> t_labels = read_binary_to_t_labels("t10k-labels.idx1-ubyte");

    discrete_test(pic_classified_in_bins,t_pictures,t_labels);
    cout << "456" << endl;
    return 0;
}