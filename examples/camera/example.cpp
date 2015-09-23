
#include <igl/Viewport.h>
#include <igl/Camera.h>
#include <igl/matlab_format.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/anttweakbar/ReAntTweakBar.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/PI.h>
#include <igl/EPS.h>
#include <igl/get_seconds.h>
#include <igl/material_colors.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/project.h>
#include <igl/opengl2/unproject.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifndef GLUT_WHEEL_UP
#define GLUT_WHEEL_UP    3
#endif
#ifndef GLUT_WHEEL_DOWN
#define GLUT_WHEEL_DOWN  4
#endif
#ifndef GLUT_WHEEL_RIGHT
#define GLUT_WHEEL_RIGHT 5
#endif
#ifndef GLUT_WHEEL_LEFT
#define GLUT_WHEEL_LEFT  6
#endif
#ifndef GLUT_ACTIVE_COMMAND
#define GLUT_ACTIVE_COMMAND 8
#endif

#include <vector>
#include <stack>
#include <iostream>
#include <algorithm>


enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type = ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP;

enum CenterType
{
  CENTER_TYPE_ORBIT = 0,
  CENTER_TYPE_FPS  = 1,
  NUM_CENTER_TYPES = 2,
} center_type = CENTER_TYPE_ORBIT;

int width,height;
#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar;
struct State
{
  std::vector<igl::Camera> cameras;
  std::vector<GLuint> tex_ids;
  std::vector<GLuint> fbo_ids;
  std::vector<GLuint> dfbo_ids;
  State():cameras(4),
    tex_ids(cameras.size()),
    fbo_ids(cameras.size()),
    dfbo_ids(cameras.size())
    {}
} s;
const Eigen::Vector4d back(1,1,1,1);
std::stack<State> undo_stack;
bool is_rotating = false;
igl::Camera down_camera;
int down_x,down_y;
std::stack<State> redo_stack;
Eigen::MatrixXd V,N;
Eigen::MatrixXi F;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

void undo()
{
  if(!undo_stack.empty())
  {
    redo_stack.push(s);
    s = undo_stack.top();
    undo_stack.pop();
  }
}

void redo()
{
  if(!redo_stack.empty())
  {
    undo_stack.push(s);
    s = redo_stack.top();
    redo_stack.pop();
  }
}

void print(const igl::Camera & camera)
{
  using namespace std;
  cout<<
    "rotation:    "<<camera.m_rotation_conj.conjugate().coeffs().transpose()<<endl<<
    "translation: "<<camera.m_translation.transpose()<<endl<<
    "eye:         "<<camera.eye().transpose()<<endl<<
    "at:          "<<camera.at().transpose()<<endl<<
    "up:          "<<camera.up().transpose()<<endl<<
    endl;
}

void init_cameras()
{
  using namespace Eigen;
  using namespace std;
  s.cameras[0].look_at(
    Vector3d(0,0,1),
    Vector3d(0,0,0),
    Vector3d(0,1,0));
  s.cameras[1].look_at(
    Vector3d(0,0,-1),
    Vector3d(0,0,0),
    Vector3d(0,1,0));
  s.cameras[2].look_at(
    Vector3d(-2,0,0),
    Vector3d(0,0,0),
    Vector3d(0,1,0));
  s.cameras[3].look_at(
    Vector3d(3,0,0),
    Vector3d(0,0,0),
    Vector3d(0,1,0));
}

bool init_render_to_texture(
  const int width, 
  const int height, 
  GLuint & tex_id, 
  GLuint & fbo_id, 
  GLuint & dfbo_id)
{
  using namespace igl;
  using namespace std;
  // Set up a "render-to-texture" frame buffer and texture combo
  glDeleteTextures(1,&tex_id);
  glDeleteFramebuffersEXT(1,&fbo_id);
  glDeleteFramebuffersEXT(1,&dfbo_id);
  // http://www.opengl.org/wiki/Framebuffer_Object_Examples#Quick_example.2C_render_to_texture_.282D.29
  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //NULL means reserve texture memory, but texels are undefined
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);
  glBindTexture(GL_TEXTURE_2D, 0);
  //-------------------------
  glGenFramebuffersEXT(1, &fbo_id);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id);
  //Attach 2D texture to this FBO
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, tex_id, 0);
  glGenRenderbuffersEXT(1, &dfbo_id);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, dfbo_id);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, width, height);
  //-------------------------
  //Attach depth buffer to FBO
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, dfbo_id);
  //-------------------------
  //Does the GPU support current FBO configuration?
  GLenum status;
  status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status)
  {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
      break;
    default:
      return false;
  }
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  return true;
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
}


// Set up double-sided lights
void lights()
{
  using namespace std;
  using namespace Eigen;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  float WHITE[4] =  {0.8,0.8,0.8,1.};
  float GREY[4] =  {0.4,0.4,0.4,1.};
  float BLACK[4] =  {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
}

void draw_scene(const igl::Camera & v_camera,
  const bool render_to_texture,
  const GLuint & v_tex_id, 
  const GLuint & v_fbo_id, 
  const GLuint & v_dfbo_id)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  if(render_to_texture)
  {
    // render to framebuffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, v_fbo_id);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, v_dfbo_id);
    glClearColor(back(0),back(1),back(2),1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMultMatrixd(v_camera.projection().data());
  //{
  //  Matrix4d m;
  //  glGetDoublev(GL_PROJECTION_MATRIX,m.data());
  //  cout<<matlab_format(m,"Camera")<<endl;
  //}
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(
    v_camera.eye()(0), v_camera.eye()(1), v_camera.eye()(2),
    v_camera.at()(0), v_camera.at()(1), v_camera.at()(2),
    v_camera.up()(0), v_camera.up()(1), v_camera.up()(2));
  //glLoadIdentity();
  //glMultMatrixd(v_camera.inverse().matrix().data());

  for(int c = 0;c<(int)s.cameras.size();c++)
  {
    auto & camera = s.cameras[c];
    if(&v_camera == &camera)
    {
      continue;
    }
    // draw camera
    glPushMatrix();
    glMultMatrixd(camera.affine().matrix().data());
    // eye
    glColor4f(0,0,0,1);
    glPointSize(10.f);
    glBegin(GL_POINTS);
    glVertex3f(0,0,0);
    glEnd();
    // frustrum
    const Vector3d u = camera.unit_plane();
    glBegin(GL_LINES);
    for(int x = -1;x<=1;x+=2)
    {
      for(int y = -1;y<=1;y+=2)
      {
        glVertex3f(0,0,0);
        glVertex3f(x*u(0),y*u(1),u(2));
      }
    }
    glEnd();
    const Vector3d n = u*(camera.m_near-FLOAT_EPS);
    glBegin(GL_QUADS);
      glVertex3f( n(0),-n(1),n(2));
      glVertex3f(-n(0),-n(1),n(2));
      glVertex3f(-n(0), n(1),n(2));
      glVertex3f( n(0), n(1),n(2));
    glEnd();
    for(int pass = 0;pass<2;pass++)
    {
      switch(pass)
      {
        case 1:
          glColor4f(1,1,1,0.5);
          glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
          //glEnable(GL_BLEND);
          glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
          glEnable(GL_TEXTURE_2D);
          glBindTexture(GL_TEXTURE_2D,s.tex_ids[c]);
          break;
        default:
        case 0:
          glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
          glColor4f(0,0,0,1);
          break;
      }
      glBegin(GL_QUADS);
        glTexCoord2d(1,0);
        glVertex3f( 0.5*u(0),-0.5*u(1),0.5*u(2));
        glTexCoord2d(0,0);
        glVertex3f(-0.5*u(0),-0.5*u(1),0.5*u(2));
        glTexCoord2d(0,1);
        glVertex3f(-0.5*u(0), 0.5*u(1),0.5*u(2));
        glTexCoord2d(1,1);
        glVertex3f( 0.5*u(0), 0.5*u(1),0.5*u(2));
      glEnd();
      switch(pass)
      {
        case 1:
          glBindTexture(GL_TEXTURE_2D, 0);
          glDisable(GL_TEXTURE_2D);
          glDisable(GL_BLEND);
          break;
        default:
          break;
      }
    }

    glPopMatrix();
  }

  // Set material properties
  lights();
  glEnable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT, GL_AMBIENT,  GOLD_AMBIENT);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,  GOLD_DIFFUSE  );
  glMaterialfv(GL_FRONT, GL_SPECULAR, GOLD_SPECULAR);
  glMaterialf (GL_FRONT, GL_SHININESS, 128);
  glMaterialfv(GL_BACK, GL_AMBIENT,  SILVER_AMBIENT);
  glMaterialfv(GL_BACK, GL_DIFFUSE,  FAST_GREEN_DIFFUSE  );
  glMaterialfv(GL_BACK, GL_SPECULAR, SILVER_SPECULAR);
  glMaterialf (GL_BACK, GL_SHININESS, 128);
  igl::opengl2::draw_mesh(V,F,N);
  glDisable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  //glLineWidth(3.f);
  //glColor4f(1,0,1,1);
  //glutWireCube(0.25);
  //glColor4f(1,0.5,0.5,1);
  ////glutWireSphere(0.125,20,20);
  {
    glPushMatrix();
    glTranslated(0,-1,0);
    igl::opengl2::draw_floor();
    glPopMatrix();
  }
  
  // Axes
  for(int d = 0;d<3;d++)
  {
    glColor4f(d==0,d==1,d==2,1);
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(d==0,d==1,d==2);
    glEnd();
  }

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  igl::opengl::report_gl_error();

  if(render_to_texture)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  }
}

void display()
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  glClearColor(back(0),back(1),back(2),0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  // Update aspect ratios (may have changed since undo/redo)
  {
    const double aspect = (double)width/(double)height;
    for(int c = 0;c<(int)s.cameras.size();c++)
    {
      auto & camera = s.cameras[c];
      auto & tex_id = s.tex_ids[c];
      auto & fbo_id = s.fbo_ids[c];
      auto & dfbo_id = s.dfbo_ids[c];
      if(aspect != camera.m_aspect)
      {
        cout<<"Initializing camera #"<<c<<"..."<<endl;
        camera.m_aspect = aspect;
        bool ret = init_render_to_texture(width,height,tex_id,fbo_id,dfbo_id);
        assert(ret);
      }
      draw_scene(camera,true,tex_id,fbo_id,dfbo_id);
    }
  }
  {
    auto & camera = s.cameras[0];
    draw_scene(camera,false,0,0,0);
  }


  TwDraw();
  glutSwapBuffers();
}


void mouse_wheel(int wheel, int direction, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  if(wheel == 0 && TwMouseMotion(mouse_x, viewport[3] - mouse_y))
  {
    static double mouse_scroll_y = 0;
    const double delta_y = 0.125*direction;
    mouse_scroll_y += delta_y;
    TwMouseWheel(mouse_scroll_y);
    return;
  }

  auto & camera = s.cameras[0];
  switch(center_type)
  {
    case CENTER_TYPE_ORBIT:
      if(wheel==0)
      {
        // factor of zoom change
        double s = (1.-0.01*direction);
        //// FOV zoom: just widen angle. This is hardly ever appropriate.
        //camera.m_angle *= s;
        //camera.m_angle = min(max(camera.m_angle,1),89);
        camera.push_away(s);
      }else
      {
        // Dolly zoom:
        camera.dolly_zoom((double)direction*1.0);
      }
      break;
    default:
    case CENTER_TYPE_FPS:
      // Move `eye` and `at` 
      camera.dolly((wheel==0?Vector3d(0,0,1):Vector3d(-1,0,0))*0.1*direction);
      break;
  }
  glutPostRedisplay();
}

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  bool tw_using = TwEventMouseButtonGLUT(glutButton,glutState,mouse_x,mouse_y);
  switch(glutButton)
  {
    case GLUT_RIGHT_BUTTON:
    case GLUT_LEFT_BUTTON:
    {
      switch(glutState)
      {
        case 1:
          // up
          glutSetCursor(GLUT_CURSOR_INHERIT);
          is_rotating = false;
          break;
        case 0:
          if(!tw_using)
          {
            push_undo();
            glutSetCursor(GLUT_CURSOR_CYCLE);
            // collect information for trackball
            is_rotating = true;
            down_camera = s.cameras[0];
            down_x = mouse_x;
            down_y = mouse_y;
          }
        break;
      }
      break;
    }
#ifdef GLUT_WHEEL_DOWN
    // Scroll down
    case GLUT_WHEEL_DOWN:
    {
      mouse_wheel(0,-1,mouse_x,mouse_y);
      break;
    }
#endif
#ifdef GLUT_WHEEL_UP
    // Scroll up
    case GLUT_WHEEL_UP:
    {
      mouse_wheel(0,1,mouse_x,mouse_y);
      break;
    }
#endif
#ifdef GLUT_WHEEL_LEFT
    // Scroll left
    case GLUT_WHEEL_LEFT:
    {
      mouse_wheel(1,-1,mouse_x,mouse_y);
      break;
    }
#endif
#ifdef GLUT_WHEEL_RIGHT
    // Scroll right
    case GLUT_WHEEL_RIGHT:
    {
      mouse_wheel(1,1,mouse_x,mouse_y);
      break;
    }
#endif
  }
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  /*bool tw_using =*/ TwMouseMotion(mouse_x,mouse_y);

  if(is_rotating)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    auto & camera = s.cameras[0];
    Quaterniond q;
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball(
          width, height,
          2.0,
          down_camera.m_rotation_conj,
          down_x, down_y, mouse_x, mouse_y,
          q);
          break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
      {
        // Rotate according to two axis valuator with fixed up vector 
        two_axis_valuator_fixed_up(
          width, height,
          2.0,
          down_camera.m_rotation_conj,
          down_x, down_y, mouse_x, mouse_y,
          q);
        break;
      }
      default:
        break;
    }
    switch(center_type)
    {
      default:
      case CENTER_TYPE_ORBIT:
        camera.orbit(q.conjugate());
        break;
      case CENTER_TYPE_FPS:
        camera.turn_eye(q.conjugate());
        break;
    }
  }
  glutPostRedisplay();
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  int mod = glutGetModifiers();
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(0);
    case 'o':
    case 'O':
      {
        s.cameras[0].m_orthographic ^= true;
        break;
      }
    case 'z':
    case 'Z':
      if(mod & GLUT_ACTIVE_COMMAND)
      {
        if(mod & GLUT_ACTIVE_SHIFT)
        {
          redo();
        }else
        {
          undo();
        }
        break;
      }
    default:
      if(!TwEventKeyboardGLUT(key,mouse_x,mouse_y))
      {
        cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
      }
  }
  glutPostRedisplay();
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;

  // print key commands
  cout<<"[Command+Z]         Undo."<<endl;
  cout<<"[Shift+Command+Z]   Redo."<<endl;
  cout<<"[^C,ESC]            Exit."<<endl;

  if(!readOFF("../shared/cheburashka.off",V,F))
  {
    cerr<<"Failed to read in mesh..."<<endl;
    return 1;
  }
  V.rowwise() -= V.colwise().minCoeff().eval();
  per_face_normals(V,F,N);

  // Init glut
  glutInit(&argc,argv);
  if( !TwInit(TW_OPENGL, NULL) )
  {
    // A fatal error occured
    fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
    return 1;
  }
  // Create a tweak bar
  rebar.TwNewBar("bar");
  TwDefine("bar label='camera' size='200 550' text=light alpha='200' color='68 68 68'");
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString(
      "RotationType","igl_trackball,two_axis_fixed_up");
  rebar.TwAddVarRW("rotation_type", RotationTypeTW,&rotation_type,
    "keyIncr=] keyDecr=[");
  TwType CenterTypeTW = igl::anttweakbar::ReTwDefineEnumFromString(
      "CenterType","orbit,fps");
  rebar.TwAddVarRW("center_type", CenterTypeTW,&center_type,
    "keyIncr={ keyDecr=}");
  rebar.TwAddVarRW("rotation", TW_TYPE_QUAT4D,s.cameras[0].m_rotation_conj.coeffs().data(),"");
  rebar.load(REBAR_NAME);
  init_cameras();

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8");
  // Top right corner
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutInitWindowPosition(glutGet(GLUT_SCREEN_WIDTH)/2.0,-1);
  glutCreateWindow("camera");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMotionFunc(mouse_drag);

  glutMainLoop();

  return 0;
}
