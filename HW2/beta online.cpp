#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

long double factorial(int x){
    if(x > 1){
        return x * factorial(x - 1);
    }else{
        return 1;
    }
}

int main()
{
    int count = 0, a = 0, b = 0;
    std::ifstream input( "testfile.txt" );
    for( std::string line; getline( input, line ); )
    {
        ////print case
        count++;
        cout << "case " << count << ": " << line << endl;
        
        ////count likelihood and print
        int n = 0,m = 0;
        for(char& c : line) {
            n++;
            if(c == '1'){
                m++;
            }
        }
        long double likelihood = (factorial(n)/(factorial(m)*factorial(n-m)))*pow(m/static_cast<long double>(n),m)*pow((n-m)/static_cast<long double>(n),n-m);
        cout << "Likelihood: " <<likelihood << endl;
        cout << "Beta prior: a = " << a << " b = " << b << endl;
        ////posterior 正比 prior * likelihood
        ////x^?            x^a   *  x^m
        ////(1-x)^?       (1-x)^b*  (1-x)^(m-n) 
        ////常係數不影響次方
        a += m;
        b += n-m;
        cout << "Beta prior: a = " << a << " b = " << b << endl << endl;
    }
    return 0;
}
