#include <igl/get_seconds.h>
#include <igl/material_colors.h>
#include <igl/png/render_to_png.h>

#ifdef __APPLE__
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

// This example displays one of the following shapes
typedef enum { SHAPE_TEAPOT=1, SHAPE_TORUS=2, SHAPE_CONE=3, SHAPE_SPHERE=4 } Shape;
#define NUM_SHAPES 4
Shape g_CurrentShape = SHAPE_TEAPOT;
int width,height;
double alpha = 0.8;
int capture_count = 0;
bool capture_on_next = false;

const float light_pos[4] = {-0.1,-0.1,1.0,0.0};

// Callback function called by GLUT to render screen
void Display(void)
{
  using namespace igl;
  using namespace std;
  float v[4]; // will be used to set light paramters
  
  // Clear frame buffer
  glClearColor(1, 1, 1, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  
  // Set light
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  v[0] = v[1] = v[2] = 0.4f; v[3] = 1.0f;
  glLightfv(GL_LIGHT0, GL_AMBIENT, v);
  v[0] = v[1] = v[2] = 0.8f; v[3] = 1.0f;
  glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
  glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
  
  // Set material
  glDisable(GL_COLOR_MATERIAL);
  float mat_ambient[4], mat_diffuse[4], mat_specular[4], mat_shininess=128;
  copy(CYAN_AMBIENT,CYAN_AMBIENT+4,mat_ambient);
  copy(CYAN_DIFFUSE,CYAN_DIFFUSE+4,mat_diffuse);
  copy(CYAN_SPECULAR,CYAN_SPECULAR+4,mat_specular);
  mat_ambient[3] =  alpha;
  mat_diffuse[3] =  alpha;
  mat_specular[3] = alpha;
  glMaterialfv(GL_BACK, GL_AMBIENT,  mat_ambient);
  glMaterialfv(GL_BACK, GL_DIFFUSE,  mat_diffuse);
  glMaterialfv(GL_BACK, GL_SPECULAR, mat_specular);
  glMaterialf( GL_BACK, GL_SHININESS, mat_shininess);

  copy(GOLD_AMBIENT,GOLD_AMBIENT+4,mat_ambient);
  copy(GOLD_DIFFUSE,GOLD_DIFFUSE+4,mat_diffuse);
  copy(GOLD_SPECULAR,GOLD_SPECULAR+4,mat_specular);
  glMaterialfv(GL_FRONT, GL_AMBIENT,  mat_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,  mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialf( GL_FRONT, GL_SHININESS, mat_shininess);
  
  if(g_CurrentShape==SHAPE_TEAPOT)
  {
    glFrontFace(GL_CW);
  }else
  {
    glFrontFace(GL_CCW);
  }
  // Rotate and draw shape
  glPushMatrix();
  glRotated(30,1,0,0);
  glRotated(360*(fmod(get_seconds(),10.0)/10.0),0,1,0);
  glScaled(1.5,1.5,1.5);
  glCallList(g_CurrentShape);
  glPopMatrix();

  if(capture_on_next)
  {
    stringstream padnum; 
    padnum << "render_to_png-example-" << setw(4) << setfill('0') << capture_count++ << ".png";
    igl::png::render_to_png(padnum.str(),width,height);
    capture_on_next = false;
  }
  
  // Present frame buffer
  glutSwapBuffers();
  
  // Recall Display at next frame
  glutPostRedisplay();
}


// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
  // Set OpenGL viewport and camera
  glViewport(0, 0, width, height);
  ::width = width;
  ::height = height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40, (double)width/height, 1, 10);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0,0,5, 0,0,0, 0,1,0);
  glTranslatef(0, 0.6f, -1);
}



void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  switch(key)
  {
    case 's':
      g_CurrentShape = (Shape)((g_CurrentShape)%NUM_SHAPES+1);
      cout<<"g_CurrentShape: "<<g_CurrentShape<<endl;
      break;
    case char(27):
      exit(0);
    case ' ':
      capture_on_next = true;
      break;
    default:
      cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
  }
}

// Main
int main(int argc, char *argv[])
{
  // Initialize GLUT
  glutInit(&argc, argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(640, 480);
  glutCreateWindow("render_to_png example (press space to render)");
  glutCreateMenu(NULL);
  
  // Set GLUT callbacks
  glutDisplayFunc(Display);
  glutReshapeFunc(Reshape);
  
  glutKeyboardFunc(key);
  
  // Create some 3D objects (stored in display lists)
  glNewList(SHAPE_TEAPOT, GL_COMPILE);
  glutSolidTeapot(1.0);
  glEndList();
  glNewList(SHAPE_TORUS, GL_COMPILE);
  glutSolidTorus(0.3, 1.0, 16, 32);
  glEndList();
  glNewList(SHAPE_CONE, GL_COMPILE);
  glutSolidCone(1.0, 1.5, 64, 4);
  glEndList();
  glNewList(SHAPE_SPHERE, GL_COMPILE);
  glutSolidSphere(1.0, 50, 40);
  glEndList();

  // Call the GLUT main loop
  glutMainLoop();
  
  return 0;
}

