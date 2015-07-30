#include <igl/C_STR.h>
#include <igl/Camera.h>
#include <igl/REDRUM.h>
#include <igl/bfs_orient.h>
#include <igl/components.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/get_seconds.h>
#include <igl/jet.h>
#include <igl/list_to_matrix.h>
#include <igl/material_colors.h>
#include <igl/normalize_row_lengths.h>
#include <igl/orient_outward.h>
#include <igl/pathinfo.h>
#include <igl/per_face_normals.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/quat_to_mat.h>
#include <igl/randperm.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/readWRL.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/unique_simplices.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/write_triangle_mesh.h>
#include <igl/anttweakbar/ReAntTweakBar.h>
#include <igl/embree/reorient_facets_raycast.h>

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

#include <ctime>
#include <string>
#include <vector>
#include <stack>
#include <iostream>

int cc_selected = -1;

Eigen::MatrixXd V;
Eigen::VectorXd Vmid,Vmin,Vmax;
double bbd = 1.0;
Eigen::MatrixXi F;
Eigen::VectorXi CC;
struct State
{
  igl::Camera camera;
  Eigen::MatrixXd N;
  Eigen::MatrixXd C;
} s;
std::string out_filename;

// See README for descriptions
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type;

enum CenterType
{
  CENTER_TYPE_ORBIT = 0,
  CENTER_TYPE_FPS  = 1,
  NUM_CENTER_TYPES = 2,
} center_type = CENTER_TYPE_ORBIT;

enum OrientMethod
{
  ORIENT_METHOD_OUTWARD = 0,
  ORIENT_METHOD_AO = 1,
  NUM_ORIENT_METHODS = 2,
} orient_method = ORIENT_METHOD_AO;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool wireframe_visible = false;
bool fill_visible = true;

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

// Forward
void init_patches();
void init_relative();

void push_undo()
{
  undo_stack.push(s);
  // Clear
  redo_stack = std::stack<State>();
}

void TW_CALL set_orient_method(const void * value, void * clientData)
{
  const OrientMethod old_orient_method = orient_method;
  orient_method = *(const OrientMethod *)value;
  if(orient_method != old_orient_method)
  {
    init_patches();
    init_relative();
  }
}

void TW_CALL get_orient_method(void * value, void *clientData)
{
  OrientMethod * om = (OrientMethod *)(value);
  *om = orient_method;
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
  glMultMatrixd(camera.projection().data());
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
  //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  float WHITE[4] = {1,1,1,1.};
  float GREY[4] = {0.4,0.4,0.4,1.};
  float BLACK[4] = {0.,0.,0.,1.};
  float NEAR_BLACK[4] =  {0.1,0.1,0.1,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,BLACK);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
  //glEnable(GL_LIGHT1);
  //pos(0) *= -1;
  //pos(1) *= -1;
  //pos(2) *= -1;
  //glLightfv(GL_LIGHT1,GL_AMBIENT,BLACK);
  //glLightfv(GL_LIGHT1,GL_DIFFUSE,NEAR_BLACK);
  //glLightfv(GL_LIGHT1,GL_SPECULAR,BLACK);
  //glLightfv(GL_LIGHT1,GL_POSITION,pos.data());
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

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();
  push_object();

  // Set material properties
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  if(wireframe_visible)
  {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    if(fill_visible)
    {
      glColor3f(0,0,0);
      igl::opengl2::draw_mesh(V,F,s.N);
    }else
    {
      igl::opengl2::draw_mesh(V,F,s.N,s.C);
    }

    // visualize selected patch
    glLineWidth(10);
    glBegin(GL_TRIANGLES);
    glColor3d(0, 0, 0);
    // loop over faces
    for(int i = 0; i<F.rows();i++)
    {
      if (CC(i) != cc_selected) continue;
      // loop over corners of triangle
      for(int j = 0;j<3;j++)
      {
        glVertex3d(V(F(i,j),0),V(F(i,j),1),V(F(i,j),2));
      }
    }
    glEnd();
    glLineWidth(1);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_TRIANGLES);
    glColor3d(1, 0, 0);
    // loop over faces
    for(int i = 0; i<F.rows();i++)
    {
      if (CC(i) != cc_selected) continue;
      // loop over corners of triangle
      glNormal3d(s.N(i,0),s.N(i,1),s.N(i,2));
      for(int j = 0;j<3;j++)
      {
        glVertex3d(V(F(i,j),0),V(F(i,j),1),V(F(i,j),2));
      }
    }
    glEnd();
  }
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  if(fill_visible)
  {
    glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
    glPolygonOffset(1.0, 0);
    igl::opengl2::draw_mesh(V,F,s.N,s.C);
  }

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
  bool tw_using = TwMouseMotion(mouse_x,mouse_y);

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
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  per_face_normals(V,F,s.N);
  normalize_row_lengths(s.N,s.N);
  Vmax = V.colwise().maxCoeff();
  Vmin = V.colwise().minCoeff();
  Vmid = 0.5*(Vmax + Vmin);
  bbd = (Vmax-Vmin).norm();
}

void randomly_color(
  const Eigen::VectorXi & CC,
  Eigen::MatrixXd & C)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  VectorXi I;
  srand ( unsigned ( time(0) ) );
  double num_cc = (double)CC.maxCoeff()+1.0;
  randperm(num_cc,I);
  C.resize(CC.rows(),3);
  for(int f = 0;f<CC.rows();f++)
  {
    jet(
      (double)I(CC(f))/num_cc,
      C(f,0),
      C(f,1),
      C(f,2));
  }
}

void TW_CALL randomize_colors(void * /*clientData*/)
{
  push_undo();
  randomly_color(CC,s.C);
}

void init_patches()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  {
    VectorXi VCC;
    components(F,VCC);
    cout<<"There are "<<VCC.maxCoeff()+1<<" connected components of vertices."<<endl;
  }
  bfs_orient(F,F,CC);
  VectorXi I;
  switch(orient_method)
  {
    case ORIENT_METHOD_AO:
    {
      cout<<"orient_outward_ao()"<<endl;
      igl::embree::reorient_facets_raycast(V,F,F,I);
      break;
    }
    case ORIENT_METHOD_OUTWARD:
    default:
      cout<<"orient_outward()"<<endl;
      orient_outward(V,F,CC,F,I);
      break;
  }
  double num_cc = (double)CC.maxCoeff()+1.0;
  cout<<"There are "<<num_cc<<" 'manifold/orientable' patches of faces."<<endl;
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

bool save(const std::string & out_filename)
{
  using namespace std;
  using namespace igl;
  if(write_triangle_mesh(out_filename,V,F))
  {
    cout<<GREENGIN("Saved mesh to `"<<out_filename<<"` successfully.")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Failed to save mesh to `"<<out_filename<<"`.")<<endl;
    return false;
  }
}

void TW_CALL saveCB(void * /*clientData*/)
{
  save(out_filename);
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  int mod = glutGetModifiers();
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(0);
    case 'I':
    case 'i':
      {
        push_undo();
        s.N *= -1.0;
        F = F.rowwise().reverse().eval();
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
      }else
      {
        push_undo();
        Quaterniond q;
        snap_to_canonical_view_quat(s.camera.m_rotation_conj,1.0,q);
        switch(center_type)
        {
          default:
          case CENTER_TYPE_ORBIT:
            s.camera.orbit(q.conjugate());
            break;
          case CENTER_TYPE_FPS:
            s.camera.turn_eye(q.conjugate());
            break;
        }
      }
      break;
    case 'u':
        mouse_wheel(0, 1,mouse_x,mouse_y);
        break;
    case 'j':
        mouse_wheel(0,-1,mouse_x,mouse_y);
        break;
    case 'n':
      cc_selected = (cc_selected + 1) % (CC.maxCoeff() + 2);
      cout << "selected cc: " << cc_selected << endl;
      glutPostRedisplay();
      break;
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
  string filename = "../shared/truck.obj";
  switch(argc)
  {
    case 3:
      out_filename = argv[2];
    case 2:
      // Read and prepare mesh
      filename = argv[1];
      break;
    default:
      cerr<<"Usage:"<<endl<<"    ./example input.obj (output.obj)"<<endl;
      cout<<endl<<"Opening default mesh..."<<endl;
      break;
  }

  // print key commands
  cout<<"[Click] and [drag]  Rotate model using trackball."<<endl;
  cout<<"[Z,z]               Snap rotation to canonical view."<<endl;
  cout<<"[Command+Z]         Undo."<<endl;
  cout<<"[Shift+Command+Z]   Redo."<<endl;
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
  }else if(ext == "ply")
  {
    // Convert extension to lower case
    if(!igl::readPLY(filename,vV,vF,vN,vTC))
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
  MatrixXi F_unique;
  unique_simplices(F, F_unique);
  F = F_unique;

  init_patches();
  init_relative();
  randomly_color(CC,s.C);

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
  TwDefine("bar label='Patches' size='200 550' text=light alpha='200' color='68 68 68'");
  rebar.TwAddVarRW("camera_rotation", TW_TYPE_QUAT4D,
    s.camera.m_rotation_conj.coeffs().data(), "open readonly=true");
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("RotationType",
    "igl_trackball,two-axis-valuator-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  TwType CenterTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("CenterType","orbit,fps");
  rebar.TwAddVarRW("center_type", CenterTypeTW,&center_type,
    "keyIncr={ keyDecr=}");
  TwType OrientMethodTW = igl::anttweakbar::ReTwDefineEnumFromString("OrientMethod",
    "outward,ambient-occlusion");
  rebar.TwAddVarCB( "orient_method", OrientMethodTW,
    set_orient_method,get_orient_method,NULL,"keyIncr=< keyDecr=>");

  rebar.TwAddVarRW("wireframe_visible",TW_TYPE_BOOLCPP,&wireframe_visible,"key=l");
  rebar.TwAddVarRW("fill_visible",TW_TYPE_BOOLCPP,&fill_visible,"key=f");
  rebar.TwAddButton("randomize_colors",randomize_colors,NULL,"key=c");
  if(out_filename != "")
  {
    rebar.TwAddButton("save",
      saveCB,NULL,
      C_STR("label='Save to `"<<out_filename<<"`' "<<
      "key=s"));
  }
  rebar.load(REBAR_NAME);


  animation_from_quat = Quaterniond(1,0,0,0);
  s.camera.m_rotation_conj = animation_from_quat;
  animation_start_time = get_seconds();

  // Init antweakbar
#ifdef __APPLE__
  glutInitDisplayString( "rgba depth double samples>=8");
#else
  glutInitDisplayString( "rgba depth double ");   // samples>=8 somehow not supported on Kenshi's machines...?
#endif
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutCreateWindow("patches");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();

  return 0;
}
