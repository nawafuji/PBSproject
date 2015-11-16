#include "ParticleObject.h"

class Balloon : public ParticleObject
{
public:
    Balloon(double r);

    Eigen::MatrixXd positions;
    double radius;


    void render();
    void update();
};
