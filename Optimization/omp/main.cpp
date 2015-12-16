#include "GLUT/glut.h"
#include <CoreServices/CoreServices.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
#include "Balloon.h"
#include "Sph.h"

#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/adjacency_list.h>


unsigned int lastTick=0;
double timestep = 0.1; //second
unsigned int stepping= (unsigned int)(timestep * 1e6); //(unsigned int)(1e6 * Scene::step);


// camera parameters
float sphi = 0.0;
float stheta = 0.0;
float sdepth = 300;
float zNear = 1.0, zFar = 1000.0;

float windowWidth = 600;
float windowHeight = 600;

int downX;
int downY;

float gridWidth = 5.0;
int gridNum = 20;

char transMode = 'e';
bool isPlaying = false;
bool once = false;

using namespace std;

Balloon* balloon = NULL;
Sph* sph = NULL;

 

unsigned int getTime()
{
    // mac
    UnsignedWide uCur;
    Microseconds(&uCur);
    return uCur.lo;
}

void init()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);

    // undisplay back side
    // glEnable(GL_CULL_FACE);
    // glCullFace(GL_BACK);

    glEnable(GL_LINE_SMOOTH);

    // Setting light
    if(false){
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      GLfloat lightpos[] = {0.0, 0.0, 50.0, 0.0};
      glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
    }

    // Setting depth
    glEnable(GL_DEPTH_TEST);

    // enable transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    lastTick = getTime();

}

void initBalloon()
{
    //obj = new Balloon(1.0);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd VN;
    std::vector<std::vector<int>> VV;

    igl::readOFF("../model/sphere_scaled.off", V, F);

    balloon = new Balloon(0.01, timestep, V, F);
}

void initSph()
{
    sph = new Sph();
    if(balloon != NULL)
    {
      sph->setRadius(balloon->getAveRadius());
      std::cout << balloon->getAveRadius() << std::endl;
    }
}

void dispGrid()
{
    glPushMatrix();
    {
      glColor3d(0.0, 0.0, 0.0);
      for(int i=-gridNum/2; i<=gridNum/2; i++){
        glLineWidth(0.1);
        glBegin(GL_LINES);
        glVertex3d(gridWidth*i, -gridWidth*gridNum/2, 0);
        glVertex3d(gridWidth*i, gridWidth*gridNum/2, 0);
        glEnd();
        glBegin(GL_LINES);
        glVertex3d(-gridWidth*gridNum/2, gridWidth*i, 0);
        glVertex3d(gridWidth*gridNum/2, gridWidth*i, 0);
        glEnd();
      }
    }
    glPopMatrix();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);//Clear the screen

    // camera transformation
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    gluPerspective(30, windowWidth/windowHeight, zNear, zFar);
    gluLookAt(0.5*sdepth,0.4*sdepth,0.3*sdepth,0,0,0,0,0,1);

    unsigned int tm = getTime();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // pseudo camera transformation
    glRotated(-stheta, 1.0, 0.0, 0.0);
    glRotated(sphi, 0.0, 0.0, 1.0);

    if (tm-lastTick > stepping)
    {
      lastTick += stepping;
      if(isPlaying)
      {
        balloon->update();
        sph->setRadius(balloon->getAveRadius());
        sph->step();
        balloon->setAirPressure(sph->getPressure()/1000000.0);
        std::cout << sph->getPressure() << std::endl;
        if(once)
        {
          isPlaying = false;
          once = false;
        }
      }
    }
    balloon->render();
    sph->render();

    dispGrid();
    
    // glFlush();//Draw everything to the screen
    glutSwapBuffers();
    glutPostRedisplay();//Start drawing again
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    // Switch matrix mode
    glMatrixMode(GL_PROJECTION);
    // Init transformation matrix
    glLoadIdentity();
    gluPerspective(30.0, (double)width / (double)height, zNear, zFar);

    // Switch matrix mode
    glMatrixMode(GL_MODELVIEW);

    windowWidth = width;
    windowHeight = height;
}

void mouse(int button, int state, int x, int y)
{
    downX = x; downY = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
    switch(transMode)
    {
      // rotate
      case 'r': 
        sphi+=(float)(x-downX)/4.0;
        stheta+=(float)(downY-y)/4.0;
        break;
      // scale
      case 'e':
        sdepth += (float)(downY - y) / 10.0;
        break;
      default:
          break;
    }
    downX = x;
    downY = y;
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key){
      // rotating mode
      case 'r':
        transMode = key;
        break;
      // scaling mode
      case 'e':
        transMode = key;
        break;
      case 'p':
        isPlaying = !isPlaying;
        break;
      case 'o':
        isPlaying = true;
        once = true;
        break;
      case 'h':
        balloon->pomp(0.02);
        break;
      case 'l':
        balloon->pomp(0.02);
        break;
      case 'b':
        balloon->burst();
        break;
        //set a new boundary sphere of size 1.5
      case 's':
        sph->setRadius(1.5);
        break;
        // double the number of particles
      case 'd':
        sph->expand();
        break;
        // print the average pressure
      case 'g':
        std::cout << sph->getPressure()<< "\n";
        break;
      default:
        break;
    }
}

void idle(void)
{
    glutPostRedisplay();
}

int main(int argc, char** argv)
{

    // Init OpenGL 
    glutInit(&argc, argv);
    // Init display mode 
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    // Init window size
    glutInitWindowSize(windowWidth, windowHeight);
    // Open the window
    glutCreateWindow("Example");
    // Set display function
    glutDisplayFunc(display);
    // Set mouse function
    glutMouseFunc(mouse);
    // Set motion function
    glutMotionFunc(motion);
    // Set keyboard function
    glutKeyboardFunc(keyboard);
    // Set reshape function
    glutReshapeFunc(reshape);
    // glutIdleFunc(idle);
    init();

    initBalloon();
    initSph();
	
    glutMainLoop();
    // delete sc;
	
	return EXIT_SUCCESS;
}
