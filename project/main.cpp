#include "EasyBMP.h"
#include "classes.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>


int main() {

    int q = 0;//q is simply an index which allows us to control whether or not the code runs.
    while (q < 1) {
        std::cout << "Do you want to investigate the chaos game(0), generalised chaos game (1) or collage theorem (2)?" << std::endl;
        int n;
        std::cin >> n;
        UI a(1920,1200);//Creates an object a which is an instance of the UI (User interface) class.
        if (n == 0) {
            a.ChaosIN();//Depending on what the user wants to investigate a different part of the user-interface class is accessed.
        }
        if (n == 1) {
            a.GChaosIN();
        }
        if (n == 2) {
            a.CollageIN();

        }
        std::cout<<"Would you like to run the program again? Yes (1) or no (0)?"<<std::endl;
        int YN;
        std::cin>>YN;
        if(YN==0){
            q++;//Program is terminated, user breaks out of the while loop
        }
        else{

        }


    }



}