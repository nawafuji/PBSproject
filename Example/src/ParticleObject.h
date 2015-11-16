#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>

class ParticleObject 
{
public:
    ParticleObject(double radius){this->radius = radius;};

    Eigen::MatrixXd positions;
    double radius;

    virtual void render(){};
    virtual void update(){};
};
