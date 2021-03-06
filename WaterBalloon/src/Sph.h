#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>
#include <random>

class Sph {
public:
    Sph();
    Eigen::MatrixX3d x;
    Eigen::MatrixX3d u;
    Eigen::MatrixXi neighbours;
    Eigen::VectorXd rho;
    Eigen::VectorXd p;
    Eigen::MatrixX3d f;
    bool useBoundary = true;
    
    int N;
    double radius_squared;

    void render();
    void update();
    void expand();
    void addParticle(int num, double v);
    void addParticles(int num, double radius);
    float getPressure();
    void setRadius(float r_new, double v);
    void discardBoundary();
    bool isFilled();

    
    void step();
    void searchNeighbour(int i);
    double w(double x_ij);
    double w_grad(double x_ij);
    double w_grad2(double x_ij);
    void computeDensity(int i);
    void computePressure(int i);
    void computeForce(int i);
};
