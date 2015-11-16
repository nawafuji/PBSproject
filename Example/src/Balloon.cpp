#include "Balloon.h"

Balloon::Balloon(double r) : ParticleObject(r)
{
    positions = Eigen::MatrixXd(3, 4);
    
    positions.col(0) << 0.0, 0.0, 5.0;
    positions.col(1) << 5.0, 5.0, 5.0;
    positions.col(2) << 5.0, 0.0, 5.0;
    positions.col(3) << 0.0, 5.0, 5.0;
    std::cout << positions << std::endl;
    radius = r;
}

void Balloon::render()
{
    glColor3f(0.5,0.5,0.5);
    for(int i=0; i<positions.cols(); i++)
    {
      glPushMatrix();
      {
        glTranslated(positions(0,i), positions(1,i), positions(2,i));
        glutSolidSphere(radius, 50, 50);
      }
      glPopMatrix();
    }
}

void Balloon::update()
{
    for(int i=0; i<positions.cols(); i++)
    {
      positions.col(i) += Eigen::Vector3d(0.01, 0.01, 0.01);
    }
}
