//TODO: reset rho(i) to 0 or rho_0 in computedensity???
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <random>
#include "Sph.h"

using  namespace Eigen;//Eigen::Matrix2Xd, Eigen::MatrixXd,Eigen::MatrixXi,Eigen::VectorXd;



const int N_start=0;
//const double d=0.05; //1.
const double d=0.02; //1.
const double h=2*d+0.00001;//2*d;
const double rho_0=1000;
const double m=d*d*d*rho_0;
const double mu=25;
const double cutOff=h;
const double boundary=0.8;

const double k=1000;

const double g=9.81;
const int N_steps=10;
const double dt=0.004;

const double radius_0=1;



Sph::Sph()
{
    N=N_start;
    radius_squared=radius_0;
    std::mt19937 mt_rand(42);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
   
    u=MatrixX3d::Zero(N,3);
    x=MatrixX3d::Zero(N,3);
    f=MatrixX3d::Zero(N,3);

    rho=VectorXd::Zero(N);
    p=VectorXd::Zero(N);
    double q=0,qq=0;
    for (int i=0; i<N; ++i) {
        x(i,0)=d*(i%10);
        if(i%10==0){++q;}
        if(i%50==0)
        {
            ++qq;
            q=0;
        }
        x(i,1)=d*(q);
        x(i,2)=d*((qq))-0.02;
        
        u(i,0)=0;//-10*dis(mt_rand);
        u(i,1)=0;//10*dis(mt_rand);
        u(i,2)=0;//10*dis(mt_rand);
        p(i)=0.;
        rho(i)=rho_0;
        
    }
     neighbours=MatrixXi::Zero(N,N);
   
}

void Sph::searchNeighbour(int i)
{
    Vector3d temp=Vector3d::Zero(3);
// #pragma omp parallel for
    for (int j=0; j< N ; ++j) { //loop over all particles
        if (j!=i)
        {
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            if(temp.norm()<cutOff){
                neighbours(i,j)=1;
                neighbours(j,i)=1;
            }else{ neighbours(i,j)=0;
                neighbours(j,i)=0;}
        }
    }
	
}



double Sph::w(double x_ij) //kernel
{
    if((x_ij<=h)&&(x_ij>=0))
    {
        return 315/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),3); //poly6
       // return 15./(M_PI*std::pow(h,6))*std::pow(std::abs(h-(x_ij)),3); //spiky
    }
    return 0.;
}

double Sph::w_grad(double x_ij) //derivative of poly6 kernel
{
  
    if((x_ij<=h)&&(x_ij>=0))
    {
        //return -315*6/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),2)*x_ij; //poly6
       // return -45./(M_PI*std::pow(h,6))*1./x_ij*std::pow((h-x_ij),2);
        return -45./(M_PI*std::pow(h,6))*std::pow((h-x_ij),2); //spiky
    }
    return 0.;
}
                       
double Sph::w_grad2(double x_ij)
{
    if((x_ij<=h)&&(x_ij>=0))
    {
       // return 315/(64*M_PI*std::pow(h,9))*(-std::pow((h*h-x_ij*x_ij),2)*6*x_ij+(h*h-x_ij*x_ij)*24*x_ij*x_ij); //poly6
       // return -90./(M_PI*std::pow(h,6))*1./(std::abs(x_ij))*(h-std::abs(x_ij))*(h-2*std::abs(x_ij)); //spiky
        return 45./(M_PI*std::pow(h,6))*(h-x_ij);
    }
    return 0.;
}

void Sph::computeDensity(int i)
{
    
    Vector3d temp=Vector3d::Zero(3);
    double rho_temp=0.;
    rho(i)=rho_0;
    for (int j=0; j<N; ++j)
    {
        if ((neighbours(i,j)==1)&& (j!=i))
        {
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            rho(i)+=m*w(temp.norm()); //check with ,neighobours' that other particle is within i.a. radius
        }

    }
   
}

void Sph::computePressure(int i)
{
    if((rho(i)-rho_0)>0.){
        p(i)=k*(rho(i)-rho_0);
    }else{p(i)=0.;}
}

void Sph::computeForce(int i)
{
    
    VectorXd f_grav(p.size());
    
    MatrixX3d f_vis=MatrixX3d::Zero(p.size(),3);
    MatrixX3d f_p=MatrixX3d::Zero(p.size(),3);
    Vector3d temp= Vector3d::Zero(3);
// #pragma omp parallel for
    for (int j=0; j<N; ++j)
    {
        if ((neighbours(i,j)==1)&& (j!=i))
        {
            //temp = Vector3d(x(i,0)-x(j,0),x(i,1)-x(j,1),x(i,2)-x(j,2));
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            
            f_p(i,1)+= m*(p(i)+p(j))/(rho(j)*2.)*(x(i,1)-x(j,1))*w_grad(temp.norm());//std::abs(x(i,1)-x(j,1)));
            f_p(i,0)+= m*(p(i)+p(j))/(rho(j)*2.)*(x(i,0)-x(j,0))*w_grad(temp.norm());
            
            f_p(i,2)+=m*(p(i)+p(j))/(rho(j)*2.)*(x(i,2)-x(j,2))*w_grad(temp.norm());
            
            
            
            f_vis(i,1)+=mu*(m/rho(j))*(u(j,1)-u(i,1))*w_grad2(temp.norm());
            f_vis(i,0)+=mu*(m/rho(j))*(u(j,0)-u(i,0))*w_grad2(temp.norm());
            f_vis(i,2)+=mu*(m/rho(j))*(u(j,2)-u(i,2))*w_grad2(temp.norm());
            
        }
        
    }
    
    //  for (int j=0;j<N; ++j) {
    f(i,1)= -f_p(i,1) +f_vis(i,1) ;//+f_grav(i);
    f(i,0)= -f_p(i,0) +f_vis(i,0); //+f_grav(i);
    f(i,2)= -f_p(i,2) +f_vis(i,2)-g*rho(i);//*rho(i);
    //}

}
                       
                       
void Sph::step()//MatrixX2d& u,MatrixX2d& x,double dt,MatrixXi& neighbours,VectorXd& p,VectorXd& rho,MatrixX2d& f)
{
    for (int i=0; i<N; ++i) {
        searchNeighbour(i);//,x,neighbours);
    }
    for (int i=0; i<N; ++i) {
        computeDensity(i);//rho,i,x,neighbours);
        computePressure(i);
    }
 
    for (int i=0; i<N; ++i)
    {
        computeForce(i);
    }
// #pragma omp parallel for
    for (int i=0; i<N; ++i) {
        u(i,1) += dt*f(i,1)/rho(i);
        u(i,0) += dt*f(i,0)/rho(i);
        u(i,2) += dt*f(i,2)/rho(i);
        double y_t =(x(i,1)+ dt*u(i,1));
        double x_t=(x(i,0)+ dt*u(i,0));
        double z_t=(x(i,2)+ dt*u(i,2));

        if(useBoundary)
        {
            // std::cout <<"sphere"<<std::endl;
          if ((x_t*x_t+y_t*y_t+z_t*z_t)>radius_squared) 
            {
             u(i,2)=-0.1*u(i,2);
               u(i,1)=-0.1*u(i,1);
              u(i,0)=-0.1*u(i,0);
            }

            else{
                x(i,1) += dt*u(i,1);
                x(i,0) += dt*u(i,0);
                x(i,2) += dt*u(i,2);

            }
        }

        if(!useBoundary){
            x(i,0) += dt*u(i,0);
            x(i,1) += dt*u(i,1);

            if((x(i,2)+ dt*u(i,2))<=-boundary) {
                x(i,2)=-boundary;
                u(i,2)=-0.1*u(i,2);
            }else if((x(i,2)+dt*u(i,2))>=boundary)
            {
                x(i,2)=boundary;
                u(i,2)=-0.1*u(i,2);
            }else
            {
                x(i,2) += dt*u(i,2);
            }
        }

    }
    
    
    
}
void Sph::render()
{
    glColor4f(0.0,0.0,0.5,0.5);
    for(int i=0; i<N; i++)
    {
        glPushMatrix();
        {
            glTranslated(x(i,0), x(i,1), x(i,2));
            glutSolidSphere(d/10, 100, 100);
        }
        glPopMatrix();
    }
}

void Sph::expand()
{
     N=N+N_steps;
    x.conservativeResize(N,3);
    u.conservativeResize(N,3);
    f.conservativeResize(N,3);
    rho.conservativeResize(N);
    neighbours.conservativeResize(N,N);
    p.conservativeResize(N);

    int q=0;
    int qq=0;
    for (int i=(N-N_steps); i<N; ++i) {
        x(i,0)=d*(i%10);
        if(i%10==0){++q;}
        if(i%50==0)
        {
            ++qq;
            q=0;
        }
        x(i,1)=d*(q);
        x(i,2)=d*((qq))+0.1;
        
        u(i,0)=0;//-10*dis(mt_rand);
        u(i,1)=0;//10*dis(mt_rand);
        u(i,2)=0;//10*dis(mt_rand);
        p(i)=0.;
        rho(i)=rho_0;
    }
}

void Sph::addParticle(int num, double v)
{
  double pressure = getPressure();
    N+=num;
    x.conservativeResize(N,3);
    u.conservativeResize(N,3);
    f.conservativeResize(N,3);
    rho.conservativeResize(N);
    neighbours.conservativeResize(N,N);
    p.conservativeResize(N);
   
    for(int i = N-num; i<N; i++){
       
        double r = (double)(rand() % 100) / 100.0 * sqrt(radius_squared);
        double theta = (double)(rand() % 100) / 100.0 * 2.0 * M_PI;

        x(i,0)=r*cos(theta);
        x(i,1)=r*sin(theta);
        x(i,2)=-sqrt(radius_squared)+0.001;
       
        u(i,0)=0;//-10*dis(mt_rand);
        u(i,1)=0;//10*dis(mt_rand);
        u(i,2)=v;//10*dis(mt_rand);
       
        rho(i)=rho_0;
        p(i) = pressure;
       
    }
   
}
void Sph::addParticles(int num, double radius)
{
  double pressure = getPressure();
    N+=num;
    x.conservativeResize(N,3);
    u.conservativeResize(N,3);
    f.conservativeResize(N,3);
    rho.conservativeResize(N);
    neighbours.conservativeResize(N,N);
    p.conservativeResize(N);
   
    for(int i = N-num; i<N; i++){
       
        double r = (double)(rand() % 100) / 100.0 * radius;
        double theta = (double)(rand() % 100) / 100.0 * 2.0 * M_PI;
        double phi = (double)(rand() % 100) / 100.0 *  M_PI;
        x(i,0)=r*sin(phi)*cos(theta);
        x(i,1)=r*sin(phi)*sin(theta);
        x(i,2)=r*cos(phi);
       
        u(i,0)=0;//-10*dis(mt_rand);
        u(i,1)=0;//10*dis(mt_rand);
        u(i,2)=0;//10*dis(mt_rand);
       
        rho(i)=rho_0;
        p(i) = pressure;
       
    }
   
}

float Sph::getPressure()
{
    float p_tot=0;
    for (int i=0;i<N;++i)
    {
        p_tot+=p(i);
    }
    return p_tot/N;
}
void Sph::setRadius(float r_new, double v)
{
    radius_squared=r_new*r_new;
    double x_t,z_t,y_t;
    for (int i=0; i<N;++i)
    {
        x_t=x(i,0);
         y_t=x(i,1);
         z_t=x(i,2);
         Eigen::Vector3d vec = x.row(i).normalized();
        if ((x_t*x_t+y_t*y_t+z_t*z_t)>radius_squared) {
            x.row(i)=vec*r_new;
            // u.row(i)=(-x.normalized()*v);
            u.row(i)=(-x.row(i).normalized()*v);
        }
    }
}

void Sph::discardBoundary(){
  setRadius(10000, 0);
  useBoundary=false;
}

bool Sph::isFilled(){

  // std::cout << "dense" << pow(d,3)*N/pow(sqrt(radius_squared),3)<<std::endl;
  return 4.0/3.0*M_PI*pow(d,3.0)*N > 4.0/3.0*M_PI*pow(sqrt(radius_squared),3.0);
}
