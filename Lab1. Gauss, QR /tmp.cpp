using namespace std;

#include <iostream>
#include <cmath>
#include <ctime>


int main (){
    srand(time(NULL));
    for (int i = 1; i <= 10; i++) 
        cout << (double)(rand() % 21 + (-10)) / 1000 << '\n';
    return 0;
}