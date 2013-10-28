// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "glutdisplay.h"
#include "sys/filename.h"
#include "lexers/streamfilters.h"
#include "lexers/parsestream.h"
#include "../tutorial_host/tutorials_host.h"

/* include GLUT for display */
#if defined(__MACOSX__)
#  include <OpenGL/gl.h>
#  include <GLUT/glut.h>
#  include <ApplicationServices/ApplicationServices.h>
#elif defined(__WIN32__)
#  include <windows.h>
#  include <GL/gl.h>   
#  include <GL/glut.h>
#else
#  include <GL/gl.h>   
#  include <GL/glut.h>
#endif

namespace embree
{
  static const double g_time0 = getSeconds();

  /* camera */
  Camera g_camera;

  /* output settings */
  static size_t g_width = 512;
  static size_t g_height = 512;

  /* fullscreen settings */
  static bool g_fullscreen = false;
  static size_t g_window_width = 512;
  static size_t g_window_height = 512;

  /* ID of created window */
  static int g_window = 0;

  /*************************************************************************************************/
  /*                                  Keyboard control                                             */
  /*************************************************************************************************/

  float g_speed = 1.0f;

  void keyboardFunc(unsigned char k, int, int)
  {
    switch (k)
    {
    case 'f' : 
      if (g_fullscreen) {
        g_fullscreen = false;
        glutReshapeWindow(g_window_width,g_window_height);
      } else {
        g_fullscreen = true;
        g_window_width = g_width;
        g_window_height = g_height;
        glutFullScreen(); 
      }
      break;
    case 'a' : g_camera.rotate(-0.02f,0.0f); break;
    case 'd' : g_camera.rotate(+0.02f,0.0f); break;
    case 'w' : g_camera.move(0.0f,0.0f,+g_speed); break;
    case 's' : g_camera.move(0.0f,0.0f,-g_speed); break;
    case '\033': case 'q': case 'Q':
      cleanup();
      glutDestroyWindow(g_window);
      exit(0);
      break;
    }
  }

  void specialFunc(int k, int, int)
  {
    switch (k) {
    case GLUT_KEY_LEFT      : g_camera.rotate(-0.02f,0.0f); break;
    case GLUT_KEY_RIGHT     : g_camera.rotate(+0.02f,0.0f); break;
    case GLUT_KEY_UP        : g_camera.move(0.0f,0.0f,+g_speed); break;
    case GLUT_KEY_DOWN      : g_camera.move(0.0f,0.0f,-g_speed); break;
    case GLUT_KEY_PAGE_UP   : g_speed *= 1.2f; break;
    case GLUT_KEY_PAGE_DOWN : g_speed /= 1.2f; break;
    }
  }

  /*************************************************************************************************/
  /*                                   Mouse control                                               */
  /*************************************************************************************************/

  static int mouseMode = 0;
  static int clickX = 0, clickY = 0;
  static bool flip14 = false;

  void clickFunc(int button, int state, int x, int y) 
  {
    if (state == GLUT_UP) {
      mouseMode = 0;
    } else {
      clickX = x; clickY = y;
      int modifiers = glutGetModifiers();
      if      (button == GLUT_LEFT_BUTTON && modifiers == GLUT_ACTIVE_CTRL) mouseMode = 1;
      else if (button == GLUT_MIDDLE_BUTTON) mouseMode = 2;
      else if (button == GLUT_RIGHT_BUTTON ) mouseMode = 3;
      else if (button == GLUT_LEFT_BUTTON  ) mouseMode = 4;

      if (flip14) {
        if (mouseMode == 4) mouseMode = 1;
        else if (mouseMode == 1) mouseMode = 4;
      }
    }
  }

  void motionFunc(int x, int y)
  {
    float dClickX = float(clickX - x), dClickY = float(clickY - y);
    clickX = x; clickY = y;

    switch (mouseMode) {
    case 1: g_camera.rotateOrbit(-0.005f*dClickX,0.005f*dClickY); break;
    case 2: break;
    case 3: g_camera.dolly(-dClickY); break;
    case 4: g_camera.rotate(-0.005f*dClickX,0.005f*dClickY); break;
    }
  }
 
  /*************************************************************************************************/
  /*                                   Window control                                              */
  /*************************************************************************************************/

  void displayFunc(void) 
  {
    AffineSpace3f pixel2world = g_camera.pixel2world(g_width,g_height);

    /* render image using ISPC */
    double t0 = getSeconds();
    int* pixels = render(g_time0-t0,
                         pixel2world.l.vx,
                         pixel2world.l.vy,
                         pixel2world.l.vz,
                         pixel2world.p);
    double dt0 = getSeconds()-t0;
    
    /* draw pixels to screen */
    glRasterPos2i(-1, 1);
    glPixelZoom(1.0f, -1.0f);
    glDrawPixels(g_width,g_height,GL_RGBA,GL_UNSIGNED_BYTE,pixels);
    glFlush();
    glutSwapBuffers();
    double dt1 = getSeconds()-t0;

    /* unmap framebuffer */
    unmap();

    /* print frame rate */
    std::ostringstream stream;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(2);
    stream << "render: ";
    stream << 1.0f/dt0 << " fps, ";
    stream << dt0*1000.0f << " ms, ";
    stream << "display: ";
    stream << 1.0f/dt1 << " fps, ";
    stream << dt1*1000.0f << " ms, ";
    stream << g_width << "x" << g_height << " pixels";
    std::cout << stream.str() << std::endl;
  }

  void reshapeFunc(int width, int height) 
  {
    resize(width,height);
    glViewport(0, 0, width, height);
    g_width = width; g_height = height;
  }
  
  void idleFunc() {
    glutPostRedisplay();
  }

  /* initialize GLUT */
  void initGlut (const std::string name, const size_t width, const size_t height, const bool fullscreen, const bool mouseMode)
  {
    g_width = width;
    g_height = height;
    resize(g_width,g_height);

    g_fullscreen = fullscreen;
    flip14 = mouseMode;
    int argc = 0; char** argv = NULL; 
    glutInit(&argc, argv);
    glutInitWindowSize((GLsizei)g_width, (GLsizei)g_height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    g_window = glutCreateWindow(name.c_str());
    if (g_fullscreen) glutFullScreen();
    glutDisplayFunc(displayFunc);
    glutIdleFunc(idleFunc);
    glutKeyboardFunc(keyboardFunc);
    glutSpecialFunc(specialFunc);
    glutMouseFunc(clickFunc);
    glutMotionFunc(motionFunc);
    glutReshapeFunc(reshapeFunc);
    glutMainLoop();
  }
}
