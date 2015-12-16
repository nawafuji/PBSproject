#include "Balloon.h"


Balloon::Balloon(double r) : ParticleObject(r)
{
    m_positions = Eigen::MatrixXd(3, 4);
    
    m_positions.col(0) << 0.0, 0.0, 5.0;
    m_positions.col(1) << 5.0, 5.0, 5.0;
    m_positions.col(2) << 5.0, 0.0, 5.0;
    m_positions.col(3) << 0.0, 5.0, 5.0;
    std::cout << m_positions << std::endl;
    m_radius = r;
}

Balloon::Balloon(double r, double dt, Eigen::MatrixXd positions, Eigen::MatrixXi faces)
    :
    ParticleObject(r)
{
    m_radius = r;
    m_step = dt;
    m_positions = positions;
    m_speeds = Eigen::MatrixXd::Zero(positions.rows(), positions.cols());
    m_forces = Eigen::MatrixXd::Zero(m_positions.rows(), m_positions.cols());
    m_num_point = m_positions.rows();
    m_mass = m_mass_all / (double)m_num_point;
    c = 2.0*sqrt(m_mass*k)/3.0;
    m_faces = faces;
    m_faceactive = Eigen::VectorXi::Zero(m_faces.rows());
    for(int i=0; i<m_faceactive.rows(); i++)
    {
        m_faceactive(i) = 1;
    }
    igl::adjacency_list(m_faces, m_adjacency_list);
    igl::triangle_triangle_adjacency(m_positions, m_faces, m_facetoface);

    std::vector<std::vector<int>> tmp_list;
    igl::vertex_triangle_adjacency(m_positions, m_faces, m_vertextoface, tmp_list);
    initEdges();
    computeNormals();


    calcAveRadius();

    m_threshold_ratio = Eigen::VectorXd::Ones(m_edges.rows()) * threshold_ratio;

    // for(int i=0; i<m_positions.rows(); i++){
    //   // std::cout << "adjacent to " << i << ": " << m_adjacency_list[i].size() << std::endl;
    // }
    // for(int i=0; i<m_positions.rows(); i++){
    //   // std::cout << "normal " << i << ": " << m_normals.norm() << std::endl;
    // }
}

void Balloon::render()
{
    glColor3f(0.5,0.5,0.5);
    for(int i=0; i<m_positions.rows(); i++)
    {
      glPushMatrix();
      {
        glTranslated(m_positions(i,0), m_positions(i,1), m_positions(i,2));
        glutSolidSphere(m_radius, 50, 50);
      }
      glPopMatrix();

      glPushMatrix();
      {
        if(false)
        {
          double nl = 5;
          glLineWidth(0.01);
          glColor3d(0.0, 1.0, 0.0);
          glBegin(GL_LINES);
          glVertex3d(m_positions(i,0), m_positions(i,1), m_positions(i,2));
          glVertex3d(m_positions(i,0) + m_normals(i,0) * nl, m_positions(i,1) + m_normals(i,1) * nl, m_positions(i,2) + m_normals(i,2) * nl);
          glEnd();
        }
      }
      glPopMatrix();
    }

    for(int i=0; i<m_faces.rows(); i++)
    {
      glPushMatrix();
      {
        if(m_faceactive(i) == 1)
        {
          glLineWidth(0.01);
          glColor4d(1.0, 0.0, 0.0, 1.0);
          if(renderTriangle)
          {
            glBegin(GL_TRIANGLES);
          }
          else{
          glBegin(GL_LINE_LOOP);
          }
          glVertex3d(m_positions(m_faces(i,0),0), m_positions(m_faces(i,0),1), m_positions(m_faces(i,0),2));
          glVertex3d(m_positions(m_faces(i,1),0), m_positions(m_faces(i,1),1), m_positions(m_faces(i,1),2));
          glVertex3d(m_positions(m_faces(i,2),0), m_positions(m_faces(i,2),1), m_positions(m_faces(i,2),2));
          glEnd();
        }
      }
      glPopMatrix();
    }
}

void Balloon::update()
{
    // update position & speed
    Eigen::MatrixXd X_n(m_positions.rows(), m_positions.cols());
    Eigen::MatrixXd V_n(m_speeds.rows(), m_speeds.cols());
    for(int i=0; i<m_positions.rows(); i++)
    {
      Eigen::Vector3d v_n;
      Eigen::Vector3d x_p = m_positions.row(i);
      Eigen::Vector3d v_p = m_speeds.row(i);
      Eigen::Vector3d f_p = m_forces.row(i);
      if(true)
      {
        // v_n = ((m_mass-m_step*(-c*m_adjacency_list[i].size()))*v_p +(f_p-c*v_p)*m_step)/(m_mass-m_step*(-c*m_adjacency_list[i].size())-m_step*m_step*(-k*m_adjacency_list[i].size()));
        // v_n = ((m_mass-m_step*(-c))*v_p + f_p)/(m_mass-m_step*(-c)-m_step*m_step*(-k));
        v_n = ((m_mass-m_step*(-d))*v_p + f_p-d*v_p)/(m_mass-m_step*(-d)-m_step*m_step*(-k));
        // v_n = ((m_mass-m_step*(-d))*v_p +(f_p-d*v_p)*m_step)/(m_mass-m_step*(-d)-m_step*m_step*(-k*m_adjacency_list[i].size()));
      }
      {
        v_n = v_p + m_step / m_mass * f_p;
      }
      X_n.row(i) = x_p + m_step * v_n;
      V_n.row(i) = v_n;
    }
    m_positions = X_n;
    m_speeds = V_n;
    checkfloor();

    calcAveRadius();

    m_forces = Eigen::MatrixXd::Zero(m_positions.rows(), m_positions.cols());

    if(!isActive)
    {
      // m_air_pressure /= 5.0;
      for(int i=0; i<m_positions.rows(); i++)
      {
        m_forces.row(i) += m_mass * g * Eigen::Vector3d(0,0,-1);
      }

      for(int i=0; i<m_faces.rows(); i++)
      {
        if(m_faceactive(i) == 0)
        {
          for(int j=0; j<3; j++)
          {
            std::vector<int> faceids = m_vertextoface[m_faces(i,j)];
            for(int k=0; k<faceids.size(); k++)
            {
              m_threshold_ratio(m_facetoedge(faceids[k],0)) /= 2.0;
              m_threshold_ratio(m_facetoedge(faceids[k],1)) /= 2.0;
              m_threshold_ratio(m_facetoedge(faceids[k],2)) /= 2.0;
            }
          }
        }
      }
    }
    // check activity
    for(int i=0; i<m_edges.rows(); i++)
    {
      if(m_edgeactive(i) == 1)
      {
        Eigen::Vector3d l = m_positions.row(m_edges(i,0)) - m_positions.row(m_edges(i,1));
        if(l.norm() > m_edgelength(i)*threshold_ratio)
        {
          m_edgeactive(i) = 0;
          m_faceactive(m_edgetoface(i,0)) = 0;
          m_faceactive(m_edgetoface(i,1)) = 0;
          isActive = false;
        }
      }
    }

    // calc forces

    for(int i=0; i<m_edges.rows(); i++)
    {
      if(m_edgeactive(i) == 1)
      {
        Eigen::Vector3d f_n = Eigen::Vector3d(0,0,0);
        Eigen::Vector3d x_1 = m_positions.row(m_edges(i,0));
        Eigen::Vector3d x_2 = m_positions.row(m_edges(i,1));
        Eigen::Vector3d v_1 = m_speeds.row(m_edges(i,0));
        Eigen::Vector3d v_2 = m_speeds.row(m_edges(i,1));

        if((x_1-x_2).norm()>m_edgelength(i)){
          f_n -= k*((x_1 - x_2).norm() - m_edgelength(i))*(x_1 - x_2).normalized();
        }
        else{
          f_n -= 1.0*k*((x_1 - x_2).norm() - m_edgelength(i))*(x_1 - x_2).normalized();
        }

        f_n -= c*((v_1 - v_2).dot((x_1 - x_2).normalized()))*(x_1 - x_2).normalized();

        m_forces.row(m_edges(i,0)) += f_n;
        m_forces.row(m_edges(i,1)) += -f_n;
      }
    }

    computeNormals();

    if(useWater)
    {
      m_forces += m_air_pressure * m_normals / (double)m_average_radius;
      // std::cout << "p_force" << m_air_pressure * (double)m_average_radius*(double)m_average_radius*4.0*M_PI/m_num_point<< std::endl;
    }
    {
      m_forces += m_air_pressure * m_normals / (double)m_average_radius;
      // std::cout << "p_force" << m_air_pressure / (double)m_average_radius << std::endl;
    }

}

void Balloon::computeNormals()
{
    Eigen::MatrixXd facenormals;
    m_normals = Eigen::MatrixXd::Zero(m_positions.rows(), m_positions.cols());
    igl::per_face_normals(m_positions, m_faces, facenormals);

    for(int i=0; i<m_positions.rows(); i++)
    {
      Eigen::Vector3d normal(0,0,0);
      for(int j=0; j<m_vertextoface[i].size(); j++)
      {
        if(m_faceactive(i) == 1)
        {
          normal += facenormals.row(m_vertextoface[i][j]);
        }
      }
      if(normal.norm() != 0)
      {
        m_normals.row(i) = normal.normalized();
      }
    }
}

void Balloon::initEdges()
{
    igl::edge_topology(m_positions, m_faces, m_edges, m_facetoedge, m_edgetoface);
    m_edgeactive = Eigen::VectorXi::Zero(m_edges.rows());
    for(int i=0; i<m_edgeactive.rows(); i++)
    {
        m_edgeactive(i) = 1;
    }
    m_edgelength = Eigen::VectorXd::Zero(m_edges.rows());
    for(int i=0; i<m_edges.rows(); i++)
    {
        Eigen::Vector3d l = m_positions.row(m_edges(i,0)) - m_positions.row(m_edges(i,1));
        m_edgelength(i) = l.norm();
    }
}

void Balloon::pomp(double p)
{
    m_air_pressure += p; 
}

/* TODO: implement */
void Balloon::burst()
{
    int id = 100;
    m_faceactive(id) = 0;
    m_edgeactive(m_facetoedge(id,0)) = 0;
    m_edgeactive(m_facetoedge(id,1)) = 0;
    m_edgeactive(m_facetoedge(id,2)) = 0;

    isActive = false;
}

void Balloon::calcAveRadius()
{
    double r = 0;
    for (int i = 0; i<m_positions.rows(); i++) {
      r+=m_positions.row(i).norm();
    }
    m_average_radius = r/(double)m_positions.rows();
}

double Balloon::getAveRadius()
{
    return m_average_radius;
}

double Balloon::getAveSpeed()
{
    return m_speeds.row(0).norm();
}

void Balloon::setAirPressure(double p)
{
    m_air_pressure = p;
}

void Balloon::pomp()
{
    double pompstep = 0.000001;
    m_air_pressure += pompstep;
}

void Balloon::checkfloor()
{
    double limit = -0.5;
    for (int i = 0; i<m_positions.rows(); i++) {
      if(m_positions(i,2) <= limit){
        m_positions(i,2) = limit;
        m_speeds.row(i) = m_speeds.row(i)/10.0;
      }
    }
}
