#include "classes.h"
#include <cstdlib>
#include "EasyBMP.h"
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
//Here is the constructor of our inputs class
UI::UI(int M1, int N1) {
    M = M1;
    N = N1;
    data = new double *[M];

    for (int i = 0; i < M; ++i) {

        data[i] = new double[N];

    }

    //Initialise pixel array to zero (i.e. empty image)
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            data[i][j] = 0;
        }
    }
}

void UI::ChaosIN(){
    std::cout << "How many verticies would you be interested in considering? Pick at least 3.  Ensure the total probability of selecting a vertex is 1, also ensure x is between 0 and 1900 and y is between 0 and 1200. Otherwise the code won't compile properly. " << std::endl;
    int Num;
    std::cin >> Num;
    std::cout << "Please enter the coordinates of vertex 1 and the associated probability of selecting it"<< std::endl;
    int x1, y1;
    double p1;//probability is a number between 0-1.
    std::cin >> x1 >> y1 >> p1;
    std::cout << "Please enter the coordinates of vertex 2 and the associated probability of selecting it"<< std::endl;
    int x2, y2;
    double p2;
    std::cin >> x2 >> y2 >> p2;
    std::cout << "Please enter the coordinates of vertex 3 and the associated probability of selecting it"<< std::endl;
    int x3, y3;
    double p3;
    std::cin >> x3 >> y3 >> p3;
    Chaos a(x1, y1, x2, y2, x3, y3, p1, p2, p3);//Create an object of type Sierpinski with the user inputs
    if (Num > 3) {//Allows the user to add more vertices if they wish
        for (int j = 0; j < (Num - 3); j++) {
            std::cout << "Please enter the coordinates of vertex " << j + 4<< " and the associated probability of selecting it" << std::endl;
            int xn, yn;
            double pn;
            std::cin >> xn >> yn >> pn;
            a.AddVertex(xn, yn, pn);
        }
    }
    std::cout<< "When a vertex is selected the algorithm moves towards this vertex by a factor of 1/n, what value of n do you want to consider? (2 is default)"<< std::endl;
    int n1;
    std::cin >> n1;
    std::cout<< "How many iterations would you like to run?"<<std::endl;
    int n2;
    std::cin >> n2;
    a.Fractals(data, n1, n2, Num);//Uses the fractals method which implements the algorithm
    OUTPUTs();//Outputs image

}

void UI::GChaosIN(){//UI for generalised Choas game
    std::cout<<"Would you like to enter the contractions to have available manually (0) or input from a file (1)?"<<std::endl;
    int MF;
    std::cin>>MF;
    if(MF==0) {
        std::cout<<"Name the file where you would like store the contractions, you can then reuse them in later iterations by passing this filename to the program."<<std::endl;
        std::string file;
        std::cin >> file;
        std::ofstream myFile;
        myFile.open(file.c_str());

        std::cout<< "add contraction, specify the value of a11,a12,a21,a22,b1,b2 and P where P is the probability of selecting the contraction (as a percentage) "<< std::endl;
        double a11, a12, a21, a22, b1, b2, P;
        std::cin >> a11 >> a12 >> a21 >> a22 >> b1 >> b2 >> P;
        myFile << a11 << " " << a12 << " " << a21 << " " << a22 << " " << b1 << " " << b2 << " " << P<< std::endl;
        int k = 0; //variable which controls whether while loop is exectuted
        while (k < 1) {
            std::cout << "Add another contraction? Yes=1, No=0" << std::endl;
            int YN;
            std::cin >> YN;
            if (YN == 0) {
                k++; //If user doesn't want to add another contraction the loop is ended
            } else if (YN == 1) {
                std::cout<< "Specify the value of a11,a12,a21,a22,b1,b2 and P where P is the probability of selecting the contraction as a percentage"<< std::endl;
                double a11n, a12n, a21n, a22n, b1n, b2n, Pn;
                std::cin >> a11n >> a12n >> a21n >> a22n >> b1n >> b2n >> Pn;
                myFile << a11n << " " << a12n << " " << a21n << " " << a22n << " " << b1n << " " << b2n << " "<< Pn << std::endl;

            }
        }
        std::cout << "How many iterations would you like to run?" << std::endl;
        int iter;
        std::cin >> iter;
        GChaos b(file.c_str());//Create an object of type chaos with the filename containing the relevant contractions passed to it.
        b.Transform(iter, data);//Now apply the transformation method associated with the chaos class
        OUTPUTs();//Call UI for outputs
    }
    else if(MF==1){
        std::cout<<"What is the name of the file of contractions you would like to input?"<<std::endl;
        std::string file;
        std::cin >> file;
        GChaos c(file);
        std::cout << "How many iterations would you like to run?" << std::endl;
        int iter;
        std::cin >> iter;
        c.Transform(iter,data);
        OUTPUTs();

    }

}
void UI::OUTPUTs(){//Once an algorithm has been implemented such that the pixels "data" have been set to certain values, this UI allows the user to control the outputted image.
    double min = 9e9;
    double max = -9e9;
    for( int i=0 ; i < M ;i++ )
    {
        for( int j=0; j <N ; j++ )
        {
            if( data[i][j] < min )
            { min = data[i][j]; }
            if( data[i][j] > max )
            { max = data[i][j]; }
        }
    }


    //Create output bmp container
    BMP Output;
    Output.SetSize(M,N);
    Output.SetBitDepth(32);

    // plot the pixels
    for( int i=0 ; i < M ; i++ )
    {
        for( int j=0; j < N ; j++ )
        {
            double scaled_value = 1 - ( data[i][j] - min )/( max-min + 1e-16 );
            ebmpBYTE pixel_value = (ebmpBYTE) ( scaled_value * 255.0 );
            Output(i,j)->Red = pixel_value;   //
            Output(i,j)->Green = pixel_value; //
            Output(i,j)->Blue = pixel_value;  //
        }
    }
    //Allows user to output image to file of their choice
    std::cout<<"where do you want to output the image? (select a file name of the form name.bmp)"<<std::endl;
    std::string file;
    std::cin>>file;
    Output.WriteToFile( file.c_str() );
    //Remove data array from heap
    for(int i = 0; i < M; ++i){

        delete[] data[i];
    }

    //Remove array container
    delete[] data;



}
void UI::CollageIN(){
    std::cout<<"Do you want to input the mappings from a file (0) or enter them manually? (1)"<<std::endl;
    int YN;
    std::cin>>YN;
    if(YN==1){
        std::ofstream myFile;
        myFile.open("Mappings.txt");
        std::cout<<"Enter your first contraction mapping, enter the parameters a11, a12, a21, a22, b1, b2"<<std::endl;
        double a_11, a_12, a_21, a_22, b_1, b_2;
        std::cin>>a_11>>a_12>>a_21>>a_22>>b_1>>b_2;
        myFile << a_11 << " " << a_12 << " " << a_21 << " " << a_22 << " " << b_1 << " " << b_2<< std::endl;
        int n=0;
        while(n<1) {
            std::cout << "Do you want to add an additional mapping? Yes(1) or No(0)." << std::endl;
            int YN;
            std::cin >> YN;
            if (YN == 1) {
                std::cout<<"Enter the parameters a11, a12, a21, a22, b1, b2"<<std::endl;
                double a_11n, a_12n, a_21n, a_22n, b_1n, b_2n;
                std::cin >> a_11n >> a_12n >> a_21n >> a_22n >> b_1n >> b_2n;
                myFile << a_11n << " " << a_12n << " " << a_21n << " " << a_22n << " " << b_1n << " " << b_2n<< std::endl;
            } else if (YN == 0) {
                n++;
            }
        }
        Collage a("Mappings.txt");
        std::cout<<"Name the image you want to experiment on"<<std::endl;
        std::string imagename;
        std::cin>>imagename;
        a.DrawCollage(imagename.c_str());
    }
    if(YN==0){
        std::cout<<"Name the file you want to use mappings from"<<std::endl;
        std::string filename;
        std::cin >> filename;
        Collage a(filename.c_str());
        std::cout<<"Name the image you want to experiment on"<<std::endl;
        std::string imagename;
        std::cin>>imagename;
        a.DrawCollage(imagename.c_str());


    }
    //Remove data array from heap
    for(int i = 0; i < M; ++i){

        delete[] data[i];
    }

    //Remove array container
    delete[] data;


}

//};

void Chaos::AddVertex(int v_x, int v_y, double p){//Allows for an additional vertex to be added to the class, with an associated probability
    vx.push_back(v_x);
    vy.push_back(v_y);
    prob.push_back(p);
}
int* Chaos::SelectVertex(int M[]){//method for selecting random vertex
    double j=rand() % 100;
    int s=0;//S is an index denoting the contraction being considered, if the random number j is greater than prob[s] then S is increased
    int i=0;
    int y=(prob[0]*100);
    while(i<1){
        if(j<=y){
            M[0]=vx[s];
            M[1]=vy[s];
            i++;
        }
        else{
            y+=(prob[s+1]*100);
            s++;
        }

    }
    return M;
}


void Chaos::Midpoint(int M[2], int n2){//method which finds next point by selecting random vertex and moving towards it by some distance
    int M1[2]={};
    int*M3=SelectVertex(M1);//vector M3 is a vertex selected using selecting vertex
    M[0]=M[0]+(M3[0]-M[0])/n2;
    M[1]=M[1]+(M3[1]-M[1])/n2;

}

void Chaos::Fractals(double** dat, int n1, int n, int Num){//This is the main method which implements the algorithm.
    double p=0;
    for(int j=0; j<Num;j++){
        p+=prob[j];
    }
    if(p<0.99||p>1){//ensure probabilities sum to 1.
        std::cout<<"Invalid probabilites, make sure they sum to 1"<<std::endl;
    }
    else {
        int i = 0;//Index to control while loop
        int M2[2] = {vx[0], vy[0]};//initial chosen vertex to start from
        while (i < n) {
            Midpoint(M2, n1);//Calls the midpoint method, gets point between current position and chosen vertex
            int j = M2[0];//x component of this point
            int k = M2[1];//y component
            dat[j - 1][k-1] = 8;//updates the pixel and this point to a nonzero value
            i++;//code is rerun

        }
    }
}


GChaos::GChaos(std::string a) {//constructor for chaos class, this takes a file containing the contractions and stores the contractions in vectors of the class
    iterator = 0; //iterator set to zero,
    std::string line;
    std::ifstream myFile(a.c_str());
    if (myFile.is_open()) {

        while (std::getline(myFile, line)) {
            std::vector<double> S;
            std::stringstream linestream(line);
            float data;
            while (linestream >> data) {
                S.push_back(data);//contractions from files added to S vector as an intermediate step

            }
            if ((pow(S[0], 2) + pow(S[2], 2)) < 1 && (pow(S[1], 2) + pow(S[3], 2)) < 1 &&
                (pow(S[0], 2) + pow(S[3], 2) + pow(S[2], 2) + pow(S[1], 2) -
                 pow((Det(S[0], S[1], S[2], S[3])), 2))<1) {//this condition ensures any contractions passed from the file are valid
                a11.push_back(S[0]);//file contents added to vectors of class via intermediate step of being transfered to S vector
                a12.push_back(S[1]);
                a21.push_back(S[2]);
                a22.push_back(S[3]);
                b1.push_back(S[4]);
                b2.push_back(S[5]);
                p.push_back(S[6]);
                iterator += 1;//iterator increased by one.
            } else {
                std::cout << "Contraction "<<iterator+1<<" is not a valid contraction" << std::endl;
                iterator +=1;
            }

        }


        myFile.close();
    }
}

double GChaos::Det(double a11,double a12,double a21,double a22){//method for finding determinant of matrix
    return (a11*a22)-(a21*a12);

}

void GChaos::Transform(int n,double** data){//This is the main method where the algorithm is implemented
    double prob=0;
    for(int j=0; j<iterator; j++){
        prob+=p[j];
    }
    if(prob<99||prob>100) {
        std::cout << "invalid probabilities, make sure they sum to 100" << std::endl;
    }
    else {
        S11 = 0;//these values essentially correspond to the maximum and minimum x and y components generated by the algorithm, this is what allows us to scale the results of the algorithm to fit the pixel display
        S12 = 0;
        S21 = 0;
        S22 = 0;
        int i = 0;//index which controls the while loop
        std::vector<double> x;//these vectors store the coordinates of the output vectors of the algorithm
        std::vector<double> y;
        double M[2] = {1, 1};//start at the point 1,1 and apply a transformation to this
        while (i < n) {
            int k = SelectContraction();//selects a contraction randomly according to probabilities
            NewVector(M, k);//applies the selected transformation to our point M to get new point M;
            x.push_back(M[0]);//adds this point to vector
            y.push_back(M[1]);
            i++;


        }
        if (S11 <= 0 && S21 <=
                        0) {//following section applies appropriate scaling to all points found by algorithm such that they fill pixel screen
            for (int i = 0; i < n; i++) {
                int a1 = (x[i] - S11) * (1900 / (S12 - S11));
                int a2 = (y[i] - S21) * (1200 / (S22 - S21));
                data[a1][a2] = 8;//points set to non-zero value, the actual value is basically arbitrary as long as it's not zero
            }
        } else if (S11 > 0 && S21 <= 0) {
            for (int i = 0; i < n; i++) {
                int a1 = (x[i]) * (1900 / S12);
                int a2 = (y[i] - S21) * (1200 / (S22 - S21));
                data[a1][a2] = 8;
            }
        } else if (S11 <= 0 && S21 > 0) {
            for (int i = 0; i < n; i++) {
                int a1 = (x[i] - S11) * (1900 / (S12 - S11));
                int a2 = (y[i]) * (1200 / S22);
                data[a1][a2] = 8;
            }
        } else {
            for (int i = 0; i < n; i++) {
                int a1 = (x[i]) * (1900 / S12);
                int a2 = (y[i]) * (1200 / S22);
                data[a1][a2] = 8;
            }

        }
        S11 = 0;
        S12 = 0;
        S21 = 0;
        S22 = 0;
    }
}

void GChaos::NewVector(double vector[], int k){//here we take our current point as an input vector and apply the selected transformation denoted by k index.
    vector[0]=(a11[k]*vector[0])+(a12[k]*vector[1])+(b1[k]);
    vector[1]=(a21[k]*vector[0])+(a22[k]*vector[1])+(b2[k]);
    if(vector[0]<S11){//Here S11,S12,S21,S22 are updated if certain conditions are met
        S11=vector[0];
    }
    if(vector[0]>S12){
        S12=vector[0];
    }
    if(vector[1]<S21){
        S21=vector[1];
    }
    if(vector[1]>S22){
        S22=vector[1];
    }

}

int GChaos::SelectContraction(){//this sections selects a contraction at random according to their probabilities
    double j=rand() % 100;
    int s=0;
    int i=0;
    int k=0;
    int y=p[0];
    while(i<1){
        if(j<y){
            k=s;
            i++;
        }
        else{
            y+=p[s+1];
            s++;
        }
    }
    return k;
}

Collage::Collage(std::string a){//Collage constructor, this takes a file containing mappings as an input and then stores them in vectors
    contractnum=0;
    std::string line;
    std::ifstream myFile(a.c_str());
    if (myFile.is_open()) {

        while (std::getline(myFile, line)) {
            std::vector<double> S;
            std::stringstream linestream(line);
            float data;
            while (linestream >> data) {
                S.push_back(data);

            }
            a11.push_back(S[0]);
            a12.push_back(S[1]);
            a21.push_back(S[2]);
            a22.push_back(S[3]);
            b1.push_back(S[4]);
            b2.push_back(S[5]);
            contractnum += 1;


        }


        myFile.close();
    }
}

void Collage::DrawCollage(std::string b){//this is the draw collage method, it takes a string containing the image the user wants to process as an input
    BMP ImageIn;
    ImageIn.ReadFromFile(b.c_str());
    const unsigned int widthPixels = ImageIn.TellWidth();
    const unsigned int heightPixels = ImageIn.TellHeight();
    std::cout<<heightPixels<<std::endl;

    double **Data = new double*[widthPixels];

    for (unsigned int xPixel = 0; xPixel < widthPixels; ++xPixel){

        Data[xPixel] = new double[heightPixels];

        for (unsigned int yPixel = 0; yPixel < heightPixels; ++yPixel){
            // convert each pixel to greyscale
            double Temp = 0.30*(ImageIn(xPixel, yPixel)->Red) + 0.59*(ImageIn(xPixel, yPixel)->Green) + 0.11*(ImageIn(xPixel, yPixel)->Blue);
            if (Temp != 255){
                Data[xPixel][yPixel] = (Temp);
            }
            else{
                Data[xPixel][yPixel] = 0;
            }
        }
    }
    for(int i=0; i<contractnum; i++) {//Here each mapping is applied to the image
        for (unsigned int xPixel = 0; xPixel < widthPixels; ++xPixel) {
            for (unsigned int yPixel = 0; yPixel < heightPixels; ++yPixel) {
                // convert each pixel to greyscale
                if (Data[xPixel][yPixel] != 0) {
                    int a = (a11[i] * xPixel) + (a21[i] * yPixel) +(b1[i]);
                    int b = (a12[i] * xPixel) + (a22[i] * yPixel) +(b2[i]);
                    Data[a][b] = (i+1)*100;//The pixels of each transformed version of the object are set to a different multiple of 100, this ensures all subleaves have different shadings and hence can be seen by the user
                } else {
                    Data[xPixel][yPixel] = 0;
                }
            }
        }
    }



    double min = 9e9;
    double max = -9e9;
    for( unsigned int i=0 ; i < widthPixels ;i++ )
    {
        for( unsigned int j=0; j <heightPixels ; j++ )
        {
            if( Data[i][j] < min )
            { min = Data[i][j]; }
            if( Data[i][j] > max )
            { max = Data[i][j]; }
        }
    }

    //Create output bmp container
    BMP Output;
    Output.SetSize(widthPixels,heightPixels);
    Output.SetBitDepth(32);

    // plot the pixels
    for( unsigned int i=0 ; i < widthPixels ; i++ )
    {
        for( unsigned int j=0; j < heightPixels ; j++ )
        {
            double scaled_value = 1 - ( Data[i][j] - min )/( max-min + 1e-16 );
            ebmpBYTE pixel_value = (ebmpBYTE) ( scaled_value * 255.0 );
            Output(i,j)->Red = pixel_value;   //
            Output(i,j)->Green = pixel_value; //for more interesting colours you can play around with these
            Output(i,j)->Blue = pixel_value;  //
        }
    }
    Output.WriteToFile( "image.bmp" );
    std::cout<<"Would you like to store the contractions as an IFS code? Yes (1) or No (0)?"<<std::endl;//allows user to store mappings as an IFS code with probabilities of their choice.
    int YN;
    std::cin>>YN;
    if(YN==1){
        std::cout<<"Name the file where you would like to store the contractions (name.txt)"<<std::endl;
        std::string filename;
        std::cin>>filename;
        std::ofstream myFile;
        myFile.open(filename.c_str());
        int i=0;
        while(i<contractnum){
            std::cout << "Give the probability of selecting contraction " << i + 1 << " as a percentage probability." << std::endl;
            double prob;
            std::cin >> prob;
            if((pow(a11[i], 2) + pow(a21[i], 2)) < 1 && (pow(a12[i], 2) + pow(a22[i], 2)) < 1 &&
               (pow(a11[i], 2) + pow(a22[i], 2) + pow(a12[i], 2) + pow(a21[i], 2) -
                pow((a11[i]*a22[i]-a12[i]*a21[i]), 2))<1){
                myFile << a11[i] << " " << a12[i] << " " << a21[i] << " " << a22[i] << " " << b1[i] << " " << b2[i] << " "<< prob<<std::endl;
                i++;
            }
            else{
                std::cout<<"contraction "<<i+1<<" is not a valid contraction"<<std::endl;
                i=contractnum;
            }


        }

    }

};








