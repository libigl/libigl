// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "sys/platform.h"
#include "sys/stl/string.h"
#include "glutdisplay.h"
#include "regression.h"

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
  /* camera settings */
  extern Vec3f g_camPos;
  extern Vec3f g_camLookAt;
  extern Vec3f g_camUp;
  extern float g_camFieldOfView;
  extern float g_camRadius;
  Handle<Device::RTCamera> createCamera(const AffineSpace3f& space);
  void clearGlobalObjects();

  /* orbit camera model */
  static AffineSpace3f g_camSpace;
  static float theta;
  static float phi;
  static float psi;
  
  /* output settings */
  extern int g_refine;  
  extern bool g_fullscreen;
  extern size_t g_width, g_height;
  extern int g_numBuffers;

  /* rendering device and global handles */
  extern Device *g_device;
  extern Handle<Device::RTRenderer> g_renderer;
  extern Handle<Device::RTToneMapper> g_tonemapper;
  extern Handle<Device::RTFrameBuffer> g_frameBuffer;
  extern Handle<Device::RTScene> g_render_scene;

  /* other stuff */
  bool g_resetAccumulation = false;

  /* regression testing */
  extern bool g_regression;

  /* ID of created window */
  static int g_window = 0;

  /*************************************************************************************************/
  /*                                  Keyboard control                                             */
  /*************************************************************************************************/

  static float g_speed = 1.0f;

  void keyboardFunc(unsigned char k, int, int)
  {
    switch (k)
    {
    case 'c' : {
      AffineSpace3f cam(g_camSpace.l,g_camSpace.p);
      std::cout << "-vp " << g_camPos.x    << " " << g_camPos.y    << " " << g_camPos.z    << " " << std::endl
                << "-vi " << g_camLookAt.x << " " << g_camLookAt.y << " " << g_camLookAt.z << " " << std::endl
                << "-vu " << g_camUp.x     << " " << g_camUp.y     << " " << g_camUp.z     << " " << std::endl;
      break;
    }
    case 'f' : glutFullScreen(); break;
    case 'r' : g_refine = !g_refine; break;
    case 't' : g_regression = !g_regression; break;
    case 'l' : g_camRadius = max(0.0f, g_camRadius-1); break;
    case 'L' : g_camRadius += 1; break;
    case '\033': case 'q': case 'Q':
      clearGlobalObjects();
      glutDestroyWindow(g_window);
      exit(0);
      break;
    }

    g_resetAccumulation = true;
  }

  void specialFunc(int k, int, int)
  {
    switch (k) {
    case GLUT_KEY_LEFT      : g_camSpace = AffineSpace3f::rotate(g_camSpace.p,g_camUp,-0.01f) * g_camSpace; break;
    case GLUT_KEY_RIGHT     : g_camSpace = AffineSpace3f::rotate(g_camSpace.p,g_camUp,+0.01f) * g_camSpace; break;
    case GLUT_KEY_UP        : g_camSpace = g_camSpace * AffineSpace3f::translate(Vec3f(0,0,g_speed)); break;
    case GLUT_KEY_DOWN      : g_camSpace = g_camSpace * AffineSpace3f::translate(Vec3f(0,0,-g_speed)); break;
    case GLUT_KEY_PAGE_UP   : g_speed *= 1.2f; std::cout << "speed = " << g_speed << std::endl; break;
    case GLUT_KEY_PAGE_DOWN : g_speed /= 1.2f; std::cout << "speed = " << g_speed << std::endl; break;
    }
    g_resetAccumulation = true;
  }

  /*************************************************************************************************/
  /*                                   Mouse control                                               */
  /*************************************************************************************************/

  static int mouseMode = 0;
  static int clickX = 0, clickY = 0;

  void clickFunc(int button, int state, int x, int y)
  {
    if (state == GLUT_UP) {
      mouseMode = 0;
      if (button == GLUT_LEFT_BUTTON && glutGetModifiers() == GLUT_ACTIVE_CTRL) {
        Handle<Device::RTCamera> camera = createCamera(AffineSpace3f(g_camSpace.l,g_camSpace.p));
        Vec3f p;
        bool hit = g_device->rtPick(camera, x / float(g_width), y / float(g_height), g_render_scene, p.x, p.y, p.z);
        if (hit) {
          Vec3f delta = p - g_camLookAt;
          Vec3f right = cross(normalize(g_camUp),normalize(g_camLookAt-g_camPos));
          Vec3f offset = dot(delta,right)*right + dot(delta,g_camUp)*g_camUp;
          g_camLookAt = p;
          g_camPos += offset;
          g_camSpace = AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp);
          g_resetAccumulation = true;
        }
      }
      else if (button == GLUT_LEFT_BUTTON && glutGetModifiers() == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT)) {
        Handle<Device::RTCamera> camera = createCamera(AffineSpace3f(g_camSpace.l,g_camSpace.p));
        Vec3f p;
        bool hit = g_device->rtPick(camera, x / float(g_width), y / float(g_height), g_render_scene, p.x, p.y, p.z);
        if (hit) {
          Vec3f v = normalize(g_camLookAt - g_camPos);
          Vec3f d = p - g_camPos;
          g_camLookAt = g_camPos + v*dot(d,v);
          g_camSpace = AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp);
          g_resetAccumulation = true;
        }
      }
    }
    else {
      if (glutGetModifiers() == GLUT_ACTIVE_CTRL) return;
      clickX = x; clickY = y;
      if      (button == GLUT_LEFT_BUTTON && glutGetModifiers() == GLUT_ACTIVE_ALT) mouseMode = 4;
      else if (button == GLUT_LEFT_BUTTON)   mouseMode = 1;
      else if (button == GLUT_MIDDLE_BUTTON) mouseMode = 2;
      else if (button == GLUT_RIGHT_BUTTON)  mouseMode = 3;
    }
  }

  void motionFunc(int x, int y)
  {
    float dClickX = float(clickX - x), dClickY = float(clickY - y);
    clickX = x; clickY = y;

    // Rotate camera around look-at point (LMB + mouse move)
    if (mouseMode == 1) {
      float angularSpeed = 0.05f / 180.0f * float(pi);
      float mapping = 1.0f;
      if (g_camUp[1] < 0) mapping = -1.0f;
      theta -= mapping * dClickX * angularSpeed;
      phi += dClickY * angularSpeed;

      if (theta < 0) theta += 2.0f * float(pi);
      if (theta > 2.0f*float(pi)) theta -= 2.0f * float(pi);
      if (phi < -1.5f*float(pi)) phi += 2.0f*float(pi);
      if (phi > 1.5f*float(pi)) phi -= 2.0f*float(pi);

      float cosPhi = cosf(phi);
      float sinPhi = sinf(phi);
      float cosTheta = cosf(theta);
      float sinTheta = sinf(theta);
      float dist = length(g_camLookAt - g_camPos);
      g_camPos = g_camLookAt + dist * Vec3f(cosPhi * sinTheta, -sinPhi, cosPhi * cosTheta);
      Vec3f viewVec = normalize(g_camLookAt - g_camPos);
      Vec3f approxUp(0.0f, 1.0f, 0.0f);
      if (phi < -0.5f*float(pi) || phi > 0.5*float(pi)) approxUp = -approxUp;
      Vec3f rightVec = normalize(cross(viewVec, approxUp));
      AffineSpace3f rotate = AffineSpace3f::rotate(viewVec, psi);
      g_camUp = xfmVector(rotate, cross(rightVec, viewVec));
    }
    // Pan camera (MMB + mouse move)
    if (mouseMode == 2) {
      float panSpeed = 0.00025f;
      float dist = length(g_camLookAt - g_camPos);
      Vec3f viewVec = normalize(g_camLookAt - g_camPos);
      Vec3f strafeVec = cross(g_camUp, viewVec);
      Vec3f deltaVec = strafeVec * panSpeed * dist * float(dClickX)
        + g_camUp * panSpeed * dist * float(-dClickY);
      g_camPos += deltaVec;
      g_camLookAt += deltaVec;
    }
    // Dolly camera (RMB + mouse move)
    if (mouseMode == 3) {
      float dollySpeed = 0.01f;
      float delta;
      if (fabsf(dClickX) > fabsf(dClickY)) delta = float(dClickX);
      else delta = float(-dClickY);
      float k = powf((1-dollySpeed), delta);
      float dist = length(g_camLookAt - g_camPos);
      Vec3f viewVec = normalize(g_camLookAt - g_camPos);
      g_camPos += dist * (1-k) * viewVec;
    }
    // Roll camera (ALT + LMB + mouse move)
    if (mouseMode == 4) {
      float angularSpeed = 0.1f / 180.0f * float(pi);
      psi -= dClickX * angularSpeed;
      Vec3f viewVec = normalize(g_camLookAt - g_camPos);
      Vec3f approxUp(0.0f, 1.0f, 0.0f);
      if (phi < -0.5f*float(pi) || phi > 0.5*float(pi)) approxUp = -approxUp;
      Vec3f rightVec = normalize(cross(viewVec, approxUp));
      AffineSpace3f rotate = AffineSpace3f::rotate(viewVec, psi);
      g_camUp = xfmVector(rotate, cross(rightVec, viewVec));
    }

    g_camSpace = AffineSpace3f::lookAtPoint(g_camPos, g_camLookAt, g_camUp);
    g_resetAccumulation = true;

  }

  /*************************************************************************************************/
  /*                                   Window control                                              */
  /*************************************************************************************************/

  double t0 = getSeconds();

  void displayFunc(void)
  {
    /* create random geometry for regression test */
    if (g_regression)
      g_render_scene = createRandomScene(g_device,1,random<int>()%1000,random<int>()%1000);

    /* set accumulation mode */
    int accumulate = g_resetAccumulation ? 0 : g_refine;
    g_resetAccumulation = false;

    /* render image */
    Handle<Device::RTCamera> camera = createCamera(AffineSpace3f(g_camSpace.l,g_camSpace.p));

    /* render into framebuffer */
    double t0 = getSeconds();
    g_device->rtRenderFrame(g_renderer,camera,g_render_scene,g_tonemapper,g_frameBuffer,accumulate);
    g_device->rtSwapBuffers(g_frameBuffer);

    /* draw image in OpenGL */
    void* ptr = g_device->rtMapFrameBuffer(g_frameBuffer);
    glRasterPos2i(-1, 1);
    glPixelZoom(1.0f, -1.0f);
    glDrawPixels((GLsizei)g_width,(GLsizei)g_height,GL_RGBA,GL_FLOAT,ptr);
    g_device->rtUnmapFrameBuffer(g_frameBuffer);
    glutSwapBuffers();
    double dt = getSeconds()-t0;

    /* measure rendering time */
    std::ostringstream stream;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(2);
    stream << 1.0f/dt << " fps, ";
    stream.precision(0);
    stream << dt*1000.0f << " ms";
    std::cout << "display " << stream.str() << std::endl;
    glutSetWindowTitle((std::string("Embree: ") + stream.str()).c_str());
  }

  void reshapeFunc(int w, int h) {
    if (g_width == size_t(w) && g_height == size_t(h)) return;
    glViewport(0, 0, w, h);
    g_width = w; g_height = h;
    g_frameBuffer = g_device->rtNewFrameBuffer("RGB_FLOAT32",w,h,g_numBuffers);
    glViewport(0, 0, (GLsizei)g_width, (GLsizei)g_height);
  }

  void idleFunc() {
    glutPostRedisplay();
  }

  void GLUTDisplay(const AffineSpace3f& camera, float s, Handle<Device::RTScene> scene)
  {
    g_camSpace = camera;
    g_speed = s;
    g_render_scene = scene;
    scene = NULL; // GLUT will never end this function, thus cleanup scene by hand

    /* initialize orbit camera model */
    Vec3f viewVec = normalize(g_camLookAt - g_camPos);
    theta = atan2f(-viewVec.x, -viewVec.z);
    phi = asinf(viewVec.y);
    Vec3f approxUp(0.0f, 1.0f, 0.0f);
    if (phi < -0.5f*float(pi) || phi > 0.5*float(pi)) approxUp = -approxUp;
    Vec3f rightVec = normalize(cross(viewVec, approxUp));
    Vec3f upUnrotated = cross(rightVec, viewVec);
    psi = atan2f(dot(rightVec, g_camUp), dot(upUnrotated, g_camUp));

    /* initialize GLUT */
    int argc = 1; char* argv = (char*)"";
    glutInit(&argc, &argv);
    glutInitWindowSize((GLsizei)g_width, (GLsizei)g_height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    g_window = glutCreateWindow("Embree");
    if (g_fullscreen) glutFullScreen();
    //glutSetCursor(GLUT_CURSOR_NONE);
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
