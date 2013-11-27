// Small GLUT application to test shadow mapping for closed shapes
//

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/readWRL.h>
#include <igl/report_gl_error.h>
#include <igl/triangulate.h>
#include <igl/readOFF.h>
#include <igl/readMESH.h>
#include <igl/draw_mesh.h>
#include <igl/draw_floor.h>
#include <igl/pathinfo.h>
#include <igl/list_to_matrix.h>
#include <igl/quat_to_mat.h>
#include <igl/per_face_normals.h>
#include <igl/material_colors.h>
#include <igl/trackball.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/REDRUM.h>
#include <igl/Camera.h>
#include <igl/ReAntTweakBar.h>
#include <igl/get_seconds.h>
#include <igl/jet.h>
#include <igl/randperm.h>
#include <igl/normalize_row_lengths.h>
#include <igl/boost/components.h>
#include <igl/boost/bfs_orient.h>
#include <igl/orient_outward.h>
#include <igl/embree/orient_outward_ao.h>
#include <igl/unique_simplices.h>
#include <igl/C_STR.h>
#include <igl/write.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef WIN32
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#ifndef GLUT_WHEEL_UP
#define GLUT_WHEEL_UP    3
#define GLUT_WHEEL_DOWN  4
#define GLUT_WHEEL_RIGHT 5
#define GLUT_WHEEL_LEFT  6
#define GLUT_ACTIVE_COMMAND 1
#endif

#include <ctime>
#include <string>
#include <vector>
#include <stack>
#include <iostream>

struct Mesh
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
};
std::vector<Mesh> meshes;

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
igl::ReTwBar rebar;

// Forward
void init_mesh();

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
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
    copy(s.camera.rotation,s.camera.rotation+4,animation_from_quat.coeffs().data());
    const Vector3d up = animation_from_quat.matrix() * Vector3d(0,1,0);
    Vector3d proj_up(0,up(1),up(2));
    if(proj_up.norm() == 0)
    {
      proj_up = Vector3d(0,1,0);
    }
    proj_up.normalize();
    Quaterniond dq;
    dq = Quaterniond::FromTwoVectors(up,proj_up);
    animation_to_quat = dq * animation_from_quat;
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
}

void push_scene()
{
  using namespace igl;
  using namespace std;
  const double angle = s.camera.angle;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  double zNear = 1e-2;
  double zFar = 100;
  double aspect = ((double)width)/((double)height);
  // Amount of scaling needed to "fix" perspective z-shift
  double z_fix = 1.0;
  // 5 is far enough to see unit "things" well
  const double camera_z = 2;
  // Test if should be using true orthographic projection
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
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(0,0,camera_z,0,0,0,0,1,0);
  // Adjust scale to correct perspective
  glScaled(z_fix,z_fix,z_fix);
  // scale, pan
  glScaled( s.camera.zoom, s.camera.zoom, s.camera.zoom);
  double mat[4*4];
  quat_to_mat(s.camera.rotation,mat);
  glMultMatrixd(mat);
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
  float WHITE[4] = {1,1,1,1.};
  float BLACK[4] = {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,BLACK);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
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
    Quaterniond q;
    q.coeffs() = 
      animation_to_quat.coeffs()*t + animation_from_quat.coeffs()*(1.-t);
    q.normalize();
    copy(q.coeffs().data(),q.coeffs().data()+4,s.camera.rotation);
  }

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();

  for(auto & mesh : meshes)
  {
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
    draw_mesh(mesh.V,mesh.F,mesh.N);
  }
  
  // Draw a nice floor
  glPushMatrix();
  {
    const float GREY[4] = {0.5,0.5,0.6,1.0};
    const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
    draw_floor(GREY,DARK_GREY);
  }
  glPopMatrix();

  pop_scene();

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
      s.camera.zoom *= (1.0+double(direction)*z_diff);
      const double min_zoom = 0.01;
      const double max_zoom = 10.0;
      s.camera.zoom = min(max_zoom,max(min_zoom,s.camera.zoom));
    }
  }else
  {
    if(!is_rotating)
    {
      // Change viewing angle (reshape will take care of adjust zoom)
      const double a_diff = 1.0;
      s.camera.angle += double(direction)*a_diff;
      const double min_angle = 15.0;
      s.camera.angle = 
        min(90.0,max(min_angle,s.camera.angle));
    }
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
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball<double>(
          width,
          height,
          2.0,
          down_camera.rotation,
          down_x,
          down_y,
          mouse_x,
          mouse_y,
          s.camera.rotation);
          break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
      {
        Quaterniond down_q;
        copy(down_camera.rotation,down_camera.rotation+4,down_q.coeffs().data());
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
        copy(q.coeffs().data(),q.coeffs().data()+4,s.camera.rotation);
        break;
      }
      default:
        break;
    }
  }
}

void init_mesh(Mesh & mesh)
{
  using namespace Eigen;
  using namespace igl;
  per_face_normals(mesh.V,mesh.F,mesh.N);
  normalize_row_lengths(mesh.N,mesh.N);
  // Rescale so bounding box fits in unit ball
  Vector3d Vmax = mesh.V.colwise().maxCoeff();
  Vector3d Vmin = mesh.V.colwise().minCoeff();
  Vector3d Vmid = 0.5*(Vmax + Vmin);
  mesh.V.rowwise() -= Vmid.transpose();
  const double bbd = (Vmax-Vmin).norm();
  mesh.V /= (bbd*0.5);
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
      }else
      {
        push_undo();
        igl::snap_to_canonical_view_quat<double>(
          s.camera.rotation,
          1.0,
          s.camera.rotation);
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
  vector<string> filenames;
  switch(argc)
  {
    case 2:
      // Read and prepare meshes
      for(int a = 1;a<argc;a++)
      {
        filenames.push_back(argv[a]);
      }
      break;
    default:
    case 1:
      cerr<<"Usage:"<<endl<<
        "    ./example input1.obj input2.obj input3.obj ..."<<endl;
      cerr<<endl<<"Opening default mesh..."<<endl;
      string filename = "../shared/truck.obj";
      filenames.push_back(filename);
      break;
  }

  // print key commands
  cout<<"[Click] and [drag]  Rotate model using trackball."<<endl;
  cout<<"[Z,z]               Snap rotation to canonical view."<<endl;
  cout<<"[Command+Z]         Undo."<<endl;
  cout<<"[Shift+Command+Z]   Redo."<<endl;
  cout<<"[^C,ESC]            Exit."<<endl;

  for(auto & filename : filenames)
  {
    meshes.push_back(Mesh());
    Mesh & mesh = meshes.back();
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
    }
    if(vV.size() > 0)
    {
      if(!list_to_matrix(vV,mesh.V))
      {
        return 1;
      }
      triangulate(vF,mesh.F);
    }
    init_mesh(mesh);
  }

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
  TwDefine("bar label='Shadow Mapping' size='200 550' text=light alpha='200' color='68 68 68'");
  rebar.TwAddVarRW("camera_rotation", TW_TYPE_QUAT4D,s.camera.rotation,"");
  TwType RotationTypeTW = ReTwDefineEnumFromString("RotationType","igl_trackball,two_axis_fixed_up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  rebar.load(REBAR_NAME);

  animation_from_quat = Quaterniond(1,0,0,0);
  copy(s.camera.rotation,s.camera.rotation+4,animation_to_quat.coeffs().data());
  animation_start_time = get_seconds();

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8");
  // Top right corner
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutInitWindowPosition(glutGet(GLUT_SCREEN_WIDTH)/2.0,-1);
  glutCreateWindow("Shadow Mapping");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();

  return 0;
}
