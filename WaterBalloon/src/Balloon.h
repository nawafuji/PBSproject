#include "ParticleObject.h"
#include <cmath>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>

class Balloon : public ParticleObject
{
public:
    Balloon(double r);
    Balloon(double r, double dt, Eigen::MatrixXd positions, Eigen::MatrixXi faces);

    bool isActive = true;
    bool useWater = false;
    bool renderTriangle = false;
    bool renderNormals = false;
    void computeNormals();
    void render();
    void update();
    void setAirPressure(double p);
    void initEdges();
    void pomp(double p);
    void burst();
    void calcAveRadius();
    double getAveRadius();
    double getAveSpeed();
    void pomp();
    void checkfloor();

private:
    Eigen::MatrixXd m_positions;
    Eigen::MatrixXi m_faces;
    Eigen::MatrixXd m_forces;
    Eigen::MatrixXd m_speeds;
    Eigen::MatrixXd m_normals;
    Eigen::MatrixXi m_edges;
    Eigen::MatrixXi m_facetoedge;
    Eigen::MatrixXi m_edgetoface;
    Eigen::MatrixXi m_facetoface;
    Eigen::VectorXd m_edgelength;
    Eigen::VectorXi m_edgeactive;
    Eigen::VectorXi m_faceactive;
    std::vector<std::vector<int>> m_adjacency_list;
    std::vector<std::vector<int>> m_vertextoface;
    double m_mass_all = 0.01;
    int m_num_point;
    double m_radius = 0.001;
    double m_average_radius;
    double m_step = 0.01;
    double m_mass;
    double k = 0.1;
    double c = 0.00;
    double d = 0.2;
    double g = 9.8;
    double threshold_ratio = 4.5;
    Eigen::VectorXd m_threshold_ratio;
    double m_air_pressure = 0.0;
};
