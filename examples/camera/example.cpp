#include <igl/Viewport.h>
#include <igl/report_gl_error.h>
#include <igl/ReAntTweakBar.h>
#include <igl/trackball.h>
#include <igl/PI.h>
#include <igl/EPS.h>
#include <igl/get_seconds.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef WIN32
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <vector>
#include <stack>
#include <iostream>

class Camera
{
  public:
    //  m_zoom   Zoom of camera lens {1}
    //  m_angle  Field of view angle in degrees {15}
    //  m_aspect  Aspect ratio {1}
    //  m_near  near clipping plane {1e-2}
    //  m_far  far clipping plane {100}
    //  m_rotation  Rotation part of rigid transformation of camera {identity}
    //  m_translation  Translation part of rigid transformation of camera
    //    {(0,0,1)}
    double m_zoom, m_angle, m_aspect, m_near, m_far;
    Eigen::Quaterniond m_rotation;
    Eigen::Vector3d m_translation;
    Camera():
      m_zoom(1), m_angle(15.0), m_aspect(1), m_near(1e-2), m_far(100),
      m_rotation(1,0,0,0),
      m_translation(0,0,1)
    {
    }

    Eigen::Vector3d eye() const
    {
      using namespace Eigen;
      Affine3d t = Affine3d::Identity();
      t.rotate(m_rotation);
      t.translate(m_translation);
      return t * Vector3d(0,0,0);
    }

    Eigen::Vector3d at() const
    {
      using namespace Eigen;
      Affine3d t = Affine3d::Identity();
      t.rotate(m_rotation);
      t.translate(m_translation);
      return t * Vector3d(0,0,-1);
    }

    Eigen::Vector3d up() const
    {
      using namespace Eigen;
      Affine3d t = Affine3d::Identity();
      t.rotate(m_rotation);
      return t * Vector3d(0,1,0);
    }

    void dolly(const double d)
    {
      using namespace Eigen;
      Vector3d dv(0,0,d);
      m_translation += m_rotation.conjugate() * dv;
    }

    void look_at(
      const Eigen::Vector3d & eye,
      const Eigen::Vector3d & at,
      const Eigen::Vector3d & up)
    {
      using namespace Eigen;
      using namespace std;
      using namespace igl;
      // http://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
      // Normalize vector from at to eye
      const Vector3d F = (eye-at).normalized();
      // Project up onto plane orthogonal to F and normalize
      const Vector3d proj_up = (up-(up.dot(F))*F).normalized();
      Quaterniond a,b;
      a.setFromTwoVectors(Vector3d(0,0,-1),-F);
      b.setFromTwoVectors(a*Vector3d(0,1,0),proj_up);
      m_rotation = a*b;
      m_translation = m_rotation.conjugate() * eye;
      assert(           (eye-this->eye()).squaredNorm() < DOUBLE_EPS);
      assert((F-(this->eye()-this->at())).squaredNorm() < DOUBLE_EPS);
      assert(        (proj_up-this->up()).squaredNorm() < DOUBLE_EPS);
    }

};

enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type = ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP;

int width,height;
#define REBAR_NAME "temp.rbr"
igl::ReTwBar rebar;
struct State
{
  int viewing_camera;
  std::vector<Camera> cameras;
  State():viewing_camera(0),cameras(2){}
} s;
std::stack<State> undo_stack;
bool is_rotating = false;
Camera down_camera;
int down_x,down_y;
std::stack<State> redo_stack;

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

void print(const Camera & camera)
{
  using namespace std;
  cout<<
    "rotation:    "<<camera.m_rotation.coeffs().transpose()<<endl<<
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
  print(s.cameras[0]);
  //s.cameras[1].look_at(
  //  Vector3d(0,0,-1),
  //  Vector3d(0,0,0),
  //  Vector3d(0,1,0));
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
}


void display()
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Update aspect ratios (may have changed since undo/redo)
  const double aspect = (double)width/(double)height;
  for(auto & camera : s.cameras)
  {
    camera.m_aspect = aspect;
  }

  //camera.m_rotation *= Quaterniond(AngleAxisd(0.01,Vector3d(0,1,0)));

  auto & camera = s.cameras[s.viewing_camera];
  const double theta = cos(2.0*PI*get_seconds()*0.1)*PI*0.05;
  const Quaterniond R(AngleAxisd(theta,Vector3d(0,1,0)));
  //// Orbit
  //camera.look_at(
  //    R*Vector3d(0,0,1),
  //    Vector3d(0,0,0),
  //    Vector3d(0,1,0));
  // First person, head rotate
  //camera.look_at(
  //    Vector3d(0,0,1),
  //    Vector3d(0,0,1)-R*Vector3d(0,0,1),
  //    Vector3d(0,1,0));

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(camera.m_angle,camera.m_aspect,camera.m_near,camera.m_far);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(
    camera.eye()(0), camera.eye()(1), camera.eye()(2),
    camera.at()(0), camera.at()(1), camera.at()(2),
    camera.up()(0), camera.up()(1), camera.up()(2));

  for(int c = 0;c<(int)s.cameras.size();c++)
  {
    // draw camera
  }

  glDisable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glLineWidth(3.f);
  glColor4f(0,0,0,1);
  glutWireCube(0.25);
  glColor4f(1,0.5,0.5,1);
  glutWireSphere(0.125,20,20);
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

  report_gl_error();
  TwDraw();
  glutSwapBuffers();
  glutPostRedisplay();
}

void mouse_wheel(int wheel, int direction, int mouse_x, int mouse_y)
{
  using namespace std;
  if(wheel == 0)
  {
    static double mouse_scroll_y = 0;
    const double delta_y = 0.125*direction;
    mouse_scroll_y += delta_y;
    // absolute scale difference when changing zooms (+1)
    const double z_diff = 0.01;
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    if(TwMouseMotion(mouse_x, viewport[3] - mouse_y))
    {
      TwMouseWheel(mouse_scroll_y);
    }else
    {
      auto & camera = s.cameras[s.viewing_camera];
      camera.dolly(double(direction)*z_diff);
      //const double min_zoom = 0.01;
      //const double max_zoom = 10.0;
      //s.camera.zoom = min(max_zoom,max(min_zoom,s.camera.zoom));
    }
  }else
  {
  }
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
            down_camera = s.cameras[s.viewing_camera];
            down_x = mouse_x;
            down_y = mouse_y;
          }
        break;
      }
      break;
      // Scroll down
      case GLUT_WHEEL_DOWN:
      {
        mouse_wheel(0,-1,mouse_x,mouse_y);
        break;
      }
      // Scroll up
      case GLUT_WHEEL_UP:
      {
        mouse_wheel(0,1,mouse_x,mouse_y);
        break;
      }
      // Scroll left
      case GLUT_WHEEL_LEFT:
      {
        mouse_wheel(1,-1,mouse_x,mouse_y);
        break;
      }
      // Scroll right
      case GLUT_WHEEL_RIGHT:
      {
        mouse_wheel(1,1,mouse_x,mouse_y);
        break;
      }
    }
  }
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
    auto & camera = s.cameras[s.viewing_camera];
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball<double>(
          width,
          height,
          2.0,
          down_camera.m_rotation.coeffs().data(),
          down_x,
          down_y,
          mouse_x,
          mouse_y,
          camera.m_rotation.coeffs().data());
          break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
      {
        Quaterniond down_q = down_camera.m_rotation;
        Vector3d axis(0,1,0);
        const double speed = 2.0;
        Quaterniond q;
        q = down_q * 
          Quaterniond(
            AngleAxisd(
              M_PI*((double)(mouse_x-down_x))/(double)width*speed/2.0,
              axis.normalized()));
        q.normalize();
        {
          Vector3d axis(1,0,0);
          const double speed = 2.0;
          if(axis.norm() != 0)
          {
            q = 
              Quaterniond(
                AngleAxisd(
                  M_PI*(mouse_y-down_y)/(double)width*speed/2.0,
                  axis.normalized())) * q;
            q.normalize();
          }
        }
        camera.m_rotation = q;
        break;
      }
      default:
        break;
    }
    const bool orbit = true;
    if(orbit)
    {
      // at should be fixed
      // Undo rotation from translation part: translation along view (from
      // `at`)
      Vector3d t = down_camera.m_rotation * down_camera.m_translation;
      // Rotate to match new rotation
      camera.m_translation = camera.m_rotation * t;
      //assert((down_camera.at() - camera.at()).squaredNorm() < DOUBLE_EPS);
    }else
    {
      // eye should be fixed
      // flip rotation?
    }
  }
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
