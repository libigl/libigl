#include <igl/get_seconds.h>
using namespace igl;
#include <cstdio>
#include <cmath>
using namespace std;

#if defined(__APPLE__)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

const int width = 1088;
const int height = 612;
int frame_counter = 0;
double start_time;
// number of frames before computing fps
const int frames_per_lap = 1000;

void Display(void)
{
  // Clear the screen with current background color
  glClearColor(
    fabs(sin(get_seconds())),
    fabs(sin(get_seconds()/3)),
    fabs(sin(get_seconds()/7)),
    0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Present frame buffer
  glutSwapBuffers();
  // Recall Display at next frame
  glutPostRedisplay();

  frame_counter++;
  if(frame_counter == frames_per_lap)
  {
    double elapsed_time = get_seconds()-start_time;
    printf("%g fps: %d frames in %g seconds\n",
      (double)frame_counter/elapsed_time,frame_counter,elapsed_time);
    // reset frame counter and timer
    start_time = get_seconds();
    frame_counter = 0;
  }
}

int main(int argc,char * argv[])
{
  // Initialize GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(width, height);
  // Center window
  glutInitWindowPosition(
      glutGet(GLUT_SCREEN_WIDTH)/2-glutGet(GLUT_INIT_WINDOW_WIDTH)/2,
      glutGet(GLUT_SCREEN_HEIGHT)/2-glutGet(GLUT_INIT_WINDOW_HEIGHT)/2);
  glutCreateWindow("");
  glutCreateMenu(NULL);
  // Set GLUT callbacks
  glutDisplayFunc(Display);
  //glutReshapeFunc(Reshape);

  // initialize timer
  start_time = get_seconds();
  // Call the GLUT main loop
  glutMainLoop();
  return 0;
}
