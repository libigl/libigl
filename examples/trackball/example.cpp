/*
  This is an example program that came with AntTweakBar modified to not use
  AntTweakBar and just show how a trackball works
  g++ -c example.cpp -o example.o 
  g++ -o example example.o -framework OpenGL -framework GLUT
  rm *.o
*/

// IGL library
#include <igl/quat_to_mat.h>
#include <igl/quat_mult.h>
#include <igl/axis_angle_to_quat.h>
#include <igl/trackball.h>
using namespace igl;

#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
using namespace std;

#if defined(_WIN32) || defined(_WIN64)
//  MiniGLUT.h is provided to avoid the need of having GLUT installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual GLUT library SDK.
#   define USE_MINI_GLUT
#endif

#if defined(USE_MINI_GLUT)
#   include "../src/MiniGLUT.h"
#elif defined(__APPLE__)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

// This example displays one of the following shapes
typedef enum { SHAPE_TEAPOT=1, SHAPE_TORUS, SHAPE_CONE } Shape;
#define NUM_SHAPES 3
Shape g_CurrentShape = SHAPE_TEAPOT;
// Shape orientation (stored as a quaternion)
float rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };
float down_rotation[4];
// Shapes material
const float g_MatAmbient[] = { 0.5f, 0.0f, 0.0f, 1.0f };
const float g_MatDiffuse[] = { 1.0f, 1.0f, 0.0f, 1.0f };
// Light parameter
const double g_LightMultiplier = 1.0f;
const float g_LightDirection[] = { -0.57735f, -0.57735f, -0.57735f };
// Keep track of mouse down
int down_mouse_x, down_mouse_y;
// keep track of size
int width, height;
float speed_factor = 1;

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  down_mouse_x = mouse_x;
  down_mouse_y = mouse_y;
  copy(rotation,rotation+4,down_rotation);
}

void mouse_move(int mouse_x, int mouse_y)
{
  trackball(
    width,
    height,
    speed_factor,
    down_rotation,
    down_mouse_x,
    down_mouse_y,
    mouse_x,
    mouse_y,
    rotation);
  glutPostRedisplay();
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  if(key == ' ')
  {
    g_CurrentShape = (Shape)((g_CurrentShape)%NUM_SHAPES+1);
  }
  glutPostRedisplay();
}

// Callback function called by GLUT to render screen
void Display(void)
{
  
  // Clear frame buffer
  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glEnable(GL_NORMALIZE);
  
  float v[4]; // will be used to set light paramters
  // Set light
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  v[0] = v[1] = v[2] = g_LightMultiplier*0.4f; v[3] = 1.0f;
  glLightfv(GL_LIGHT0, GL_AMBIENT, v);
  v[0] = v[1] = v[2] = g_LightMultiplier*0.8f; v[3] = 1.0f;
  glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
  v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
  glLightfv(GL_LIGHT0, GL_POSITION, v);
  
  // Set material
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, g_MatAmbient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, g_MatDiffuse);
  
  // Rotate and draw shape
  glPushMatrix();
  float mat[4*4]; // rotation matrix
  quat_to_mat(rotation,mat);
  glMultMatrixf(mat);
  glCallList(g_CurrentShape);
  glPopMatrix();
  
  // Present frame buffer
  glutSwapBuffers();
  
  //// Recall Display at next frame
  //glutPostRedisplay();
}


// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
  // keep track of size
  ::width = width;
  ::height = height;
  // Set OpenGL viewport and camera
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40, (double)width/height, 1, 10);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0,0,5, 0,0,0, 0,1,0);
}


// Function called at exit
void Terminate(void)
{ 
  glDeleteLists(SHAPE_TEAPOT, NUM_SHAPES);
}

// Main
int main(int argc, char *argv[])
{
  bool help_and_quit = false;
  if(argc > 1)
  {
    if(
      argv[1][0] == '-' &&
      argv[1][1] == 'h')
    {
      help_and_quit = true;
    }else
    {
      int count = sscanf(argv[1],"%g",&speed_factor);
      if(count != 1)
      {
        printf("Error: %s is not a valid speed factor.",
          argv[1]);
        help_and_quit = true;
      }
    }
  }
  if(help_and_quit)
  {
    printf(
      "Usage:\n  ./example\nor\n  ./example [positive speed factor]\n\n"
      "Interaction:\n  Drag on screen to use trackball\n"
      "  Press SPACE to change shape\n");
    return 1;
  }
  
  // Initialize GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(640, 480);
  glutCreateWindow("Simple GLUT example with trackball");
  glutCreateMenu(NULL);
  
  // Set GLUT callbacks
  glutDisplayFunc(Display);
  glutReshapeFunc(Reshape);
  atexit(Terminate);  // Called after glutMainLoop ends
  
  // Set GLUT event callbacks
  // Directly redirect GLUT mouse button events to trackball example
  glutMouseFunc(mouse);
  // Directly redirect GLUT mouse motion events to trackball example
  glutMotionFunc(mouse_move);
  //glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  // Catch keyboard action
  glutKeyboardFunc(key);
  //glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
  //TwGLUTModifiersFunc(glutGetModifiers);
  
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
  
  // Init rotation
  float axis[] = { 0.7f, 0.7f, 0.0f }; 
  float angle = 0.8f;
  axis_angle_to_quat(axis,angle,rotation);

  // Call the GLUT main loop
  glutMainLoop();
  
  return 0;
}

