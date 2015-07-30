#include <igl/Camera.h>
#include <igl/REDRUM.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/get_seconds.h>
#include <igl/list_to_matrix.h>
#include <igl/material_colors.h>
#include <igl/pathinfo.h>
#include <igl/per_face_normals.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/quat_to_mat.h>
#include <igl/random_points_on_mesh.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readWRL.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/anttweakbar/ReAntTweakBar.h>

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

#include <string>
#include <vector>
#include <stack>
#include <iostream>


Eigen::MatrixXd V,N,S;
Eigen::VectorXd Vmid,Vmin,Vmax;
double bbd = 1.0;
Eigen::MatrixXi F;
struct State
{
  igl::Camera camera;
} s;

// See README for descriptions
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool is_rotating = false;
int down_x,down_y;
igl::Camera down_camera;

bool is_animating = false;
double animation_start_time = 0;
double ANIMATION_DURATION = 0.5;
Eigen::Quaterniond animation_from_quat;
Eigen::Quaterniond animation_to_quat;

int width,height;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar;

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

// No-op setter, does nothing
void TW_CALL no_op(const void * /*value*/, void * /*clientData*/)
{
}

void TW_CALL set_rotation_type(const void * value, void * clientData)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  const RotationType old_rotation_type = rotation_type;
  rotation_type = *(const RotationType *)(value);
  if(rotation_type == ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP &&
    old_rotation_type != ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
  {
    push_undo();
    animation_from_quat = s.camera.m_rotation_conj;
    snap_to_fixed_up(animation_from_quat,animation_to_quat);
    // start animation
    animation_start_time = get_seconds();
    is_animating = true;
  }
}
void TW_CALL get_rotation_type(void * value, void *clientData)
{
  RotationType * rt = (RotationType *)(value);
  *rt = rotation_type;
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
  s.camera.m_aspect = (double)width/(double)height;
}

void push_scene()
{
  using namespace igl;
  using namespace std;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  auto & camera = s.camera;
  gluPerspective(camera.m_angle,camera.m_aspect,camera.m_near,camera.m_far);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(
    camera.eye()(0), camera.eye()(1), camera.eye()(2),
    camera.at()(0), camera.at()(1), camera.at()(2),
    camera.up()(0), camera.up()(1), camera.up()(2));
}

void push_object()
{
  using namespace igl;
  glPushMatrix();
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-Vmid(0),-Vmid(1),-Vmid(2));
}

void pop_object()
{
  glPopMatrix();
}

void pop_scene()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  using namespace Eigen;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  float WHITE[4] =  {0.8,0.8,0.8,1.};
  float GREY[4] =  {0.4,0.4,0.4,1.};
  float BLACK[4] =  {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
  pos(0) *= -1;
  pos(1) *= -1;
  pos(2) *= -1;
  glLightfv(GL_LIGHT1,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT1,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT1,GL_POSITION,pos.data());
}

void display()
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(is_animating)
  {
    double t = (get_seconds() - animation_start_time)/ANIMATION_DURATION;
    if(t > 1)
    {
      t = 1;
      is_animating = false;
    }
    Quaterniond q = animation_from_quat.slerp(t,animation_to_quat).normalized();
    auto & camera = s.camera;
    camera.orbit(q.conjugate());
  }

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();
  push_object();

  // Set material properties
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

  glPointSize(3.);
  glEnable(GL_POINT_SMOOTH);
  glColor4f(0,0.5,1,1);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for(int s = 0;s<S.rows();s++)
  {
    glVertex3d(S(s,0),S(s,1),S(s,2));
  }
  glEnd();

  pop_object();

  // Draw a nice floor
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-Vmid(1));
  glTranslated(0,floor_offset,0);
  const float GREY[4] = {0.5,0.5,0.6,1.0};
  const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
  igl::opengl2::draw_floor(GREY,DARK_GREY);
  glPopMatrix();

  pop_scene();

  igl::opengl::report_gl_error();

  TwDraw();
  glutSwapBuffers();
  glutPostRedisplay();
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
  push_undo();

  auto & camera = s.camera;
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
            down_camera = s.camera;
            down_x = mouse_x;
            down_y = mouse_y;
          }
        break;
      }
      break;
    }
    // Scroll down
    case 3:
    {
      mouse_wheel(0,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll up
    case 4:
    {
      mouse_wheel(0,1,mouse_x,mouse_y);
      break;
    }
    // Scroll left
    case 5:
    {
      mouse_wheel(1,-1,mouse_x,mouse_y);
      break;
    }
    // Scroll right
    case 6:
    {
      mouse_wheel(1,1,mouse_x,mouse_y);
      break;
    }
  }
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  if(is_rotating)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    Quaterniond q;
    auto & camera = s.camera;
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball<double>(
          width,
          height,
          2.0,
          down_camera.m_rotation_conj.coeffs().data(),
          down_x,
          down_y,
          mouse_x,
          mouse_y,
          q.coeffs().data());
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
    camera.orbit(q.conjugate());
  }
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  per_face_normals(V,F,N);
  Vmax = V.colwise().maxCoeff();
  Vmin = V.colwise().minCoeff();
  Vmid = 0.5*(Vmax + Vmin);
  bbd = (Vmax-Vmin).norm();
}

void init_samples()
{
  using namespace Eigen;
  using namespace igl;
  const int n = V.rows()*10;
  SparseMatrix<double> B;
  VectorXi FI;
  random_points_on_mesh(n,V,F,B,FI);
  S = B*V;
}

void undo()
{
  using namespace std;
  if(!undo_stack.empty())
  {
    redo_stack.push(s);
    s = undo_stack.top();
    undo_stack.pop();
  }
}

void redo()
{
  using namespace std;
  if(!redo_stack.empty())
  {
    undo_stack.push(s);
    s = redo_stack.top();
    redo_stack.pop();
  }
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  int mod = glutGetModifiers();
  const bool command_down = GLUT_ACTIVE_COMMAND & mod;
  const bool shift_down = GLUT_ACTIVE_SHIFT & mod;
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
      if(command_down)
      {
        if(shift_down)
        {
          redo();
        }else
        {
          undo();
        }
        break;
      }else
      {
        push_undo();
        Quaterniond q;
        snap_to_canonical_view_quat(s.camera.m_rotation_conj,1.0,q);
        s.camera.orbit(q.conjugate());
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
  string filename = "../shared/cheburashka.off";
  if(argc < 2)
  {
    cerr<<"Usage:"<<endl<<"    ./example input.obj"<<endl;
    cout<<endl<<"Opening default mesh..."<<endl;
  }else
  {
    // Read and prepare mesh
    filename = argv[1];
  }

  // print key commands
  cout<<"[Click] and [drag]  Rotate model using trackball."<<endl;
  cout<<"[Z,z]               Snap rotation to canonical view."<<endl;
  cout<<"[⌘ Z]               Undo."<<endl;
  cout<<"[⇧ ⌘ Z]             Redo."<<endl;
  cout<<"[^C,ESC]            Exit."<<endl;

  // dirname, basename, extension and filename
  string d,b,ext,f;
  pathinfo(filename,d,b,ext,f);
  // Convert extension to lower case
  transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  vector<vector<double > > vV,vN,vTC;
  vector<vector<int > > vF,vFTC,vFN;
  if(ext == "obj")
  {
    // Convert extension to lower case
    if(!igl::readOBJ(filename,vV,vTC,vN,vF,vFTC,vFN))
    {
      return 1;
    }
  }else if(ext == "off")
  {
    // Convert extension to lower case
    if(!igl::readOFF(filename,vV,vF,vN))
    {
      return 1;
    }
  }else if(ext == "wrl")
  {
    // Convert extension to lower case
    if(!igl::readWRL(filename,vV,vF))
    {
      return 1;
    }
  //}else
  //{
  //  // Convert extension to lower case
  //  MatrixXi T;
  //  if(!igl::readMESH(filename,V,T,F))
  //  {
  //    return 1;
  //  }
  //  //if(F.size() > T.size() || F.size() == 0)
  //  {
  //    boundary_facets(T,F);
  //  }
  }
  if(vV.size() > 0)
  {
    if(!list_to_matrix(vV,V))
    {
      return 1;
    }
    polygon_mesh_to_triangle_mesh(vF,F);
  }

  init_relative();
  init_samples();

  // Init glut
  glutInit(&argc,argv);
  if( !TwInit(TW_OPENGL, NULL) )
  {
    // A fatal error occured
    fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
    return 1;
  }
  // Create a tweak bar
  rebar.TwNewBar("TweakBar");
  rebar.TwAddVarRW("camera_rotation", TW_TYPE_QUAT4D,
    s.camera.m_rotation_conj.coeffs().data(), "open readonly=true");
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("RotationType",
    "igl_trackball,two-a...-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  rebar.load(REBAR_NAME);

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutCreateWindow("upright");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();

  return 0;
}
