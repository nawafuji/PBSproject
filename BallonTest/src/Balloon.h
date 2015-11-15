#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>

class Balloon
{
public:
    Balloon(double r);

    Eigen::MatrixXd positions;
    double radius;


    void render();
    void update();
};
