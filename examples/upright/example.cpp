// Small GLUT application to help quickly orient a model in an upright
// position, scaled and shifted to fit the unit sphere.
//
// Demonstrates trackball mouse camera control and undo/redo stack with
// immutable state object (i.e. not with reverse operations).
//
// Reads (V,F) from .obj,.off,.wrl files and saves (V,F) to .obj,.off. Any
// normals, texture coordinates, etc. in the input file will be ignored and
// lost in the output file.
//
#include <igl/read_triangle_mesh.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/pathinfo.h>
#include <igl/list_to_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/material_colors.h>
#include <igl/trackball.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/write_triangle_mesh.h>
#include <igl/REDRUM.h>

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


struct State
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  Eigen::VectorXd Vmax,Vmin,Vmid;
  double bbd;
  Eigen::Quaterniond rot;
} s;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool trackball_on = false;
int down_x,down_y;
Eigen::Quaterniond down_rot;

int width,height;
std::string out_filename;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

void reshape(int width, int height)
{
  ::width = width;
  ::height = height;
  glViewport(0,0,width,height);
}

void push_scene()
{
  using namespace igl;
  using namespace std;
  const double angle = 15;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double zNear = 1e-2;
  double zFar = 100;
  double aspect = ((double)width)/((double)height);
  // Amount of scaling needed to "fix" perspective z-shift
  double z_fix = 1.0;
  // 5 is far enough to see unit "things" well
  const double camera_z = 2;
  // Test if should be using true orthographic igl::opengl2::projection
  if(angle == 0)
  {
    glOrtho(
      -0.5*camera_z*aspect,
      0.5*camera_z*aspect,
      -0.5*camera_z,
      0.5*camera_z,
      zNear,
      zFar);
  }else
  {
    // Make sure aspect is sane
    aspect = aspect < 0.01 ? 0.01 : aspect;
    gluPerspective(angle,aspect,zNear,zFar);
    z_fix = 2.*tan(angle/2./360.*2.*M_PI);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0,0,camera_z,0,0,0,0,1,0);
  // Adjust scale to correct perspective
  glScaled(z_fix,z_fix,z_fix);
  // scale, pan
}

void push_object()
{
  glPushMatrix();
  glMultMatrixd(Eigen::Affine3d(s.rot).matrix().data());
  glScaled(2./s.bbd,2./s.bbd,2./s.bbd);
  glTranslated(-s.Vmid(0),-s.Vmid(1),-s.Vmid(2));
}

void pop_object()
{
  glPopMatrix();
}

void pop_scene()
{
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
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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


  igl::opengl2::draw_mesh(s.V,s.F,s.N);

  glDisable(GL_LIGHTING);
  glPushMatrix();
  glLineWidth(1.0);
  glColor3f(0,0,0);
  glTranslated( s.Vmid(0), s.Vmid(1), s.Vmid(2));
  glScaled(
    (s.Vmax(0)-s.Vmin(0)),
    (s.Vmax(1)-s.Vmin(1)),
    (s.Vmax(2)-s.Vmin(2)));
  glutWireCube(1.0);
  glPopMatrix();
  pop_object();

  glDisable(GL_LIGHTING);
  glPushMatrix();
  glLineWidth(2.0);
  glColor3f(1,0,0);
  //glTranslated(s.Vmid(0), s.Vmid(1), s.Vmid(2));
  glScaled(
    2*(s.Vmax(0)-s.Vmin(0))/s.bbd,
    2*(s.Vmax(1)-s.Vmin(1))/s.bbd,
    2*(s.Vmax(2)-s.Vmin(2))/s.bbd);
  glutWireCube(1.0);
  glPopMatrix();

  //glPushMatrix();
  //glRotated(90,1,0,0);
  //glColor3f(0.5,0.5,0.5);
  //glutWireSphere(1.0,20,10);
  //glPopMatrix();

  glEnable(GL_LIGHTING);
  glPushMatrix();
  glTranslated(0,-1,0);
  igl::opengl2::draw_floor();
  glPopMatrix();
  pop_scene();

  glutSwapBuffers();
}

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  switch(glutState)
  {
    case 1:
      // up
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
      trackball_on = false;
      break;
    case 0:
      push_undo();
      glutSetCursor(GLUT_CURSOR_CYCLE);
      // collect information for trackball
      trackball_on = true;
      down_rot = s.rot;
      down_x = mouse_x;
      down_y = mouse_y;
    break;
  }
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;

  if(trackball_on)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    // Rotate according to trackball
    trackball<double>(
      width,
      height,
      2.0,
      down_rot.coeffs().data(),
      down_x,
      down_y,
      mouse_x,
      mouse_y,
      s.rot.coeffs().data());
  }
  glutPostRedisplay();
}

// Bake rotation and scale into V
void bake()
{
  s.V.col(0).array() -= s.Vmid(0);
  s.V.col(1).array() -= s.Vmid(1);
  s.V.col(2).array() -= s.Vmid(2);
  s.V *= 2./s.bbd;
  s.V = (s.V * s.rot.matrix().transpose()).eval();
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  per_face_normals(s.V,s.F,s.N);
  s.rot = Quaterniond(1,0,0,0);
  s.Vmax = s.V.colwise().maxCoeff();
  s.Vmin = s.V.colwise().minCoeff();
  s.Vmid = 0.5*(s.Vmax + s.Vmin);
  s.bbd = (s.Vmax-s.Vmin).norm();
}


bool save()
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  if(write_triangle_mesh(out_filename,s.V,s.F))
  {
    cout<<GREENGIN("Current baked model written to "+out_filename+".")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Current baked model failed to write to "+out_filename+".")<<endl;
    return false;
  }
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
  const int mod = glutGetModifiers();
  const bool command_down = GLUT_ACTIVE_COMMAND & mod;
  const bool shift_down = GLUT_ACTIVE_SHIFT & mod;
  switch(key)
  {
    // Ctrl-c and esc exit
    case char(3):
    case char(27):
      exit(0);
    case 'B':
    case 'b':
      push_undo();
      bake();
      init_relative();
      break;
    case 'I':
    case 'i':
      push_undo();
      s.F = s.F.rowwise().reverse().eval();
      break;
    case 'R':
    case 'r':
      push_undo();
      init_relative();
      break;
    case 'S':
    case 's':
      save();
      break;
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
        igl::snap_to_canonical_view_quat<double>(
          s.rot.coeffs().data(),
          1.0,
          s.rot.coeffs().data());
        break;
      }
    default:
      cout<<"Unknown key command: '"<<key<<"' ("<<int(key)<<")"<<endl;
  }

  glutPostRedisplay();
}

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  if(argc < 3)
  {
    cerr<<"Usage:"<<endl<<"    ./upright input.obj output.obj"<<endl;
    return 1;
  }

  // print key commands
  cout<<"[Click] and [drag]  Rotate model using trackball."<<endl;
  cout<<"[Z,z]               Snap rotation to canonical view."<<endl;
  cout<<"[B,b]               \"Bake\" current scale, shift, and rotation "
    "into model."<<endl;
  cout<<"[⌘ Z]               Undo."<<endl;
  cout<<"[⇧ ⌘ Z]             Redo."<<endl;
  cout<<"[R,r]               Reset rotation."<<endl;
  cout<<"[S,s]               Save to output model path."<<endl;
  cout<<"[I,i]               Flip orientation of faces (gold front, green "
    "back)."<<endl;
  cout<<"[^C,ESC]            Exit."<<endl;

  // Read and prepare mesh
  string filename = argv[1];
  out_filename = argv[2];
  if(!read_triangle_mesh(filename,s.V,s.F))
  {
    cout<<"Failed to read from "<<filename<<"."<<endl;
    return EXIT_FAILURE;
  }
  init_relative();

  // Init glut
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("upright");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  //glutPassiveMotionFunc(mouse_move);
  glutMainLoop();

  return EXIT_SUCCESS;
}
