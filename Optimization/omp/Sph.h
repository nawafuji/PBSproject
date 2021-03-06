#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>
#include <random>

class Sph {
public:
    Sph();
  //  int pppp;
    Eigen::MatrixX3d x;
    Eigen::MatrixX3d u;
    Eigen::MatrixXi neighbours;
    Eigen::VectorXd rho;
    Eigen::VectorXd p;
    Eigen::MatrixX3d f;
    
    int N;
    double radius_squared;

    void render();
    void update();
    void expand();
    float getPressure();
    void setRadius(float r_new);

    
    void step();
    void searchNeighbour(int i);
    double w(double x_ij);
    double w_grad(double x_ij);
    double w_grad2(double x_ij);
    void computeDensity(int i);
    void computePressure(int i);
    void computeForce(int i);
    
    
};
