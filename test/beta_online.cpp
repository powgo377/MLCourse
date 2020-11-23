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
    // for (int j = 1;j<47;j++){
        int x = 1;
        for (int i = 1;i<47;i++){
            // int count = 1;
            // int mod = i;
            // while(mod != 1){
            //     mod = mod*i;
            //     mod = mod%47;
            //     count ++;
            // }
            x*=13;
            x = x % 47;
            cout << i<<":"<<x<<endl;

            // int c = 55;
            //     cout <<101%23<< endl;
            
        }
    // }
    return 0;
}
