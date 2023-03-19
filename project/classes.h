
#include "EasyBMP.h"
#include <vector>
#include <algorithm>
#ifndef PROJECT_CLASSES_H
#define PROJECT_CLASSES_H

class UI{//User interface class, this contains several methods which are themselves a user interface which allows the user to investigate a particular part of the script
public:
    UI(int N, int M);
    void ChaosIN();
    void GChaosIN();
    void CollageIN();
    void OUTPUTs();
private:
    double **data;//Pixel data stored
    int N;
    int M;

};

class Chaos{//This class allows the user to investigate the Sirepenski triangle and related shapes
public:
    Chaos(int v_1x,int v_1y,int v_2x,int v_2y, int v_3x, int v_3y, double p1, double p2, double p3) {
        vx.push_back(v_1x);
        vx.push_back(v_2x);
        vx.push_back(v_3x);
        vy.push_back(v_1y);
        vy.push_back(v_2y);
        vy.push_back(v_3y);
        prob.push_back(p1);
        prob.push_back(p2);
        prob.push_back(p3);
    }//constructor allows the user to add the coordinates of 3 verticies each with an associated probability
    void AddVertex(int v_x, int v_y, double prob);//Allows user to add additional vertex
    int* SelectVertex(int M[]);//Method which selects vertex at random according to probabilities
    void Midpoint(int M[2], int n2);//Method finds midpoint between current position and randomly selected vertex
    void Fractals(double** dat, int n1, int n, int Num);

protected:
    std::vector<int> vx;//vector containing all x components of verticies
    std::vector<int> vy;//vector containing all y components
    std::vector<double> prob;//Associated probability of each vertex.



};


class GChaos{//class for generalised chaos game
public:
    GChaos(std::string a);//constructor
    void Transform(int n, double** data);//method which implements the algorithm, all remaining methods are called to from this method
    void NewVector(double vector[],int k);
    double Det(double a11,double a12,double a21, double a22);
    int SelectContraction();
protected:
    std::vector<double> a11;//a11-p are all vectors containing each component of the IFS codes.
    std::vector<double> a12;
    std::vector<double> a21;
    std::vector<double> a22;
    std::vector<double> b1;
    std::vector<double> b2;
    std::vector<double> p;
    double S11;//These are just variables used to scale the results of our algorithm to fit the full pixel display
    double S12;
    double S21;
    double S22;
    int iterator;//variable tracks how many contractions have been added to the class.



};

class Collage{//Collage theorem class
public:
    Collage(std::string a); //constructor, adds mappings to the class
    void DrawCollage(std::string b); //draws the collage
protected:
    int contractnum; //number of mappings addded
    std::vector<double> a11; //mapping components stored in vectors
    std::vector<double> a12;
    std::vector<double> a21;
    std::vector<double> a22;
    std::vector<double> b1;
    std::vector<double> b2;



};


#endif //PROJECT_CLASSES_H
