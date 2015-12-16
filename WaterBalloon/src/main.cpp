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
double timestep = 0.004; //second
unsigned int stepping= (unsigned int)(timestep * 1e6); //(unsigned int)(1e6 * Scene::step);

// switch useWater to change mode
// bool useWater = false;
bool useWater = true;

// camera parameters
float sphi = 0.0;
float stheta = 0.0;
float sdepth = 3.0;
float zNear = 0.0001, zFar = 100.0;

float windowWidth = 600;
float windowHeight = 600;

int downX;
int downY;

float gridWidth = 0.05;
int gridNum = 20;

char transMode = 'e';
bool isPlaying = true;
bool isPomping = false;
bool once = false;

int wcount = 0;
double limit = -0.8;
bool renderTriangle = false;


Balloon* balloon = NULL;
Sph* sph = NULL;


using namespace std;

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

    igl::readOFF("../model/sphere_m.off", V, F);
    // igl::readOFF("../model/bunny.off", V, F);
    // igl::readOFF("../model/bumpy_cube.off", V, F);
    // igl::readOFF("../model/cat.off", V, F);

    balloon = new Balloon(0.001, timestep, V, F);
    balloon->useWater = useWater;
}

void initSph()
{
    sph = new Sph();
    if(balloon != NULL)
    {
      sph->setRadius(balloon->getAveRadius(), 0);
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
        glVertex3d(gridWidth*i, -gridWidth*gridNum/2, limit);
        glVertex3d(gridWidth*i, gridWidth*gridNum/2, limit);
        glEnd();
        glBegin(GL_LINES);
        glVertex3d(-gridWidth*gridNum/2, gridWidth*i, limit);
        glVertex3d(gridWidth*gridNum/2, gridWidth*i, limit);
        glEnd();
      }
    }
    glPopMatrix();
}

void render_string(float x, float y, std::string const& str){
    float z = 0.5f;
    //glRasterPos3f(x, y, z);
    glWindowPos2i(20, windowHeight-40);
    for(int i=0; i<str.length(); i++){
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
    }
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

    std::cout << "ddd" << std::endl;
    if (tm-lastTick > stepping)
    {
      lastTick += stepping;
      if(isPlaying)
      {
        if(isPomping)
        {
          if(useWater)
          {
            if(balloon->isActive)
            {
              if(wcount == 1)
              {
                sph->addParticle(1, 5.0);
                // sph->addParticles(1, balloon->getAveRadius());
                wcount = 0;
                // render_string(0.0f, 0.0f, "Pumping");
              }
              wcount++;
            }
            else
            {
              sph->discardBoundary();
            }
          }
          else
          {
            balloon->pomp();
            // render_string(0.0f, 0.0f, "Pumping");
          }
        }
        if(useWater)
        {
          sph->step();
          double p = sph->getPressure();
            if(sph->isFilled()){
              balloon->setAirPressure((p-1e5)/1e10/2);
            }
            else{
              balloon->setAirPressure(0);
            }
        }
        balloon->update();

        if(balloon->isActive)
        {
          sph->setRadius(balloon->getAveRadius(), balloon->getAveSpeed());
        }
        if(once)
        {
          isPlaying = false;
          once = false;
        }
      }
    }
    balloon->render();
    if(useWater)
    {
      sph->render();
    }

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
        isPomping = !isPomping;
        break;
      case 'l':
        balloon->pomp(0.02);
        break;
      case 'b':
        balloon->burst();
        break;
        // set a new boundary sphere of size 1.5
      case 's':
        // sph->setRadius(1.5);
        break;
        // double the number of particles
      case 'd':
        sph->addParticles(10, balloon->getAveRadius());
        break;
        // print the average pressure
      case 't':
        renderTriangle = !renderTriangle;
        balloon->renderTriangle = renderTriangle;
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
