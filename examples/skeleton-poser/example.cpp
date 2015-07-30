#include <igl/Camera.h>
#include <igl/opengl2/MouseController.h>
#include <igl/REDRUM.h>
#include <igl/STR.h>
#include <igl/barycenter.h>
#include <igl/bone_parents.h>
#include <igl/boundary_conditions.h>
#include <igl/boundary_facets.h>
#include <igl/centroid.h>
#include <igl/colon.h>
#include <igl/opengl2/draw_beach_ball.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/opengl2/draw_skeleton_3d.h>
#include <igl/opengl2/draw_skeleton_vector_graphics.h>
#include <igl/forward_kinematics.h>
#include <igl/get_seconds.h>
#include <igl/lbs_matrix.h>
#include <igl/material_colors.h>
#include <igl/next_filename.h>
#include <igl/normalize_row_sums.h>
#include <igl/pathinfo.h>
#include <igl/per_face_normals.h>
#include <igl/quat_to_mat.h>
#include <igl/readDMAT.h>
#include <igl/readTGF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/volume.h>
#include <igl/winding_number.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/writeTGF.h>
#include <igl/anttweakbar/ReAntTweakBar.h>
#include <igl/bbw/bbw.h>
#include <igl/cgal/remesh_self_intersections.h>
#include <igl/tetgen/mesh_with_skeleton.h>
#include <igl/tetgen/tetrahedralize.h>

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
#include <queue>
#include <stack>
#include <iostream>
#include <iomanip>

#define VERBOSE

enum SkelStyleType
{
  SKEL_STYLE_TYPE_3D = 0,
  SKEL_STYLE_TYPE_VECTOR_GRAPHICS = 1,
  NUM_SKEL_STYLE_TYPE = 2
}skel_style;

Eigen::MatrixXd V,N,W,M;
Eigen::Vector3d Vmid;
double bbd = 1.0;
Eigen::MatrixXi F;
igl::Camera camera;

Eigen::MatrixXd C;
Eigen::MatrixXi BE;
Eigen::VectorXi P,RP;

struct State
{
  igl::opengl2::MouseController mouse;
  Eigen::MatrixXf colors;
} s;

bool wireframe = false;

// See README for descriptions
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type = ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool is_rotating = false;
bool centroid_is_visible = true;
int down_x,down_y;
igl::Camera down_camera;
std::string output_weights_filename,output_pose_prefix;

struct CameraAnimation
{
  bool is_animating = false;
  double DURATION = 0.5;
  double start_time = 0;
  Eigen::Quaterniond from_quat,to_quat;
} canim;

typedef std::vector<
Eigen::Quaterniond,
  Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;

struct PoseAnimation
{
  bool is_animating = false;
  double DURATION = 2;
  double start_time = 0;
  RotationList pose;
} panim;

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
    canim.from_quat = camera.m_rotation_conj;
    snap_to_fixed_up(canim.from_quat,canim.to_quat);
    // start animation
    canim.start_time = get_seconds();
    canim.is_animating = true;
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
  camera.m_aspect = (double)width/(double)height;
  s.mouse.reshape(width,height);
}

void push_scene()
{
  using namespace igl;
  using namespace std;
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
  const float back[4] = {0.75, 0.75, 0.75,0};
  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(canim.is_animating)
  {
    double t = (get_seconds() - canim.start_time)/canim.DURATION;
    if(t > 1)
    {
      t = 1;
      canim.is_animating = false;
    }
    Quaterniond q = canim.from_quat.slerp(t,canim.to_quat).normalized();
    camera.orbit(q.conjugate());
  }

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();

  // Draw a nice floor
  glEnable(GL_DEPTH_TEST);
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-Vmid(1));
  glTranslated(0,floor_offset,0);
  const float GREY[4] = {0.5,0.5,0.6,1.0};
  const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  igl::opengl2::draw_floor(GREY,DARK_GREY);
  glDisable(GL_CULL_FACE);
  glPopMatrix();

  push_object();

  const auto & draw_skeleton = [](const MatrixXd & T)
  {
    switch(skel_style)
    {
      default:
      case SKEL_STYLE_TYPE_3D:
      {
        igl::opengl2::draw_skeleton_3d(C,BE,T,s.colors,bbd*0.5);
        break;
      }
      case SKEL_STYLE_TYPE_VECTOR_GRAPHICS:
        igl::opengl2::draw_skeleton_vector_graphics(C,BE,T);
        break;
    }
  };
  // Set material properties
  glDisable(GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT, GL_AMBIENT,GOLD_AMBIENT);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,GOLD_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR,GOLD_SPECULAR);
  glMaterialf (GL_FRONT, GL_SHININESS, 128);
  glMaterialfv(GL_BACK, GL_AMBIENT,SILVER_AMBIENT);
  glMaterialfv(GL_BACK, GL_DIFFUSE,FAST_GREEN_DIFFUSE);
  glMaterialfv(GL_BACK, GL_SPECULAR,SILVER_SPECULAR);
  glMaterialf (GL_BACK, GL_SHININESS, 128);
  if(wireframe)
  {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  }
  glLineWidth(1.0);
  MatrixXd T;
  RotationList dQ;
  if(panim.is_animating)
  {
    double t = (get_seconds() - panim.start_time)/panim.DURATION;
    if(t > 1)
    {
      t = 1;
      panim.is_animating = false;
    }
    const auto & ease = [](const double t)
    {
      return 3.*t*t-2.*t*t*t;
    };
    double f = (t<0.5?ease(2.*t):ease(2.-2.*t));
    dQ.resize(panim.pose.size());
    for(int e = 0;e<(int)panim.pose.size();e++)
    {
      dQ[e] = panim.pose[e].slerp(f,Quaterniond::Identity()).normalized();
    }
  }else
  {
    dQ = s.mouse.rotations();
  }
  forward_kinematics(C,BE,P,dQ,T);
  MatrixXd U = M*T;
  MatrixXd UN;
  per_face_normals(U,F,UN);
  igl::opengl2::draw_mesh(U,F,UN);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  glDisable(GL_DEPTH_TEST);
  draw_skeleton(T);

  if(centroid_is_visible)
  {
    Vector3d cen;
    centroid(U,F,cen);
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    glTranslated(cen(0),cen(1),cen(2));
    glScaled(bbd/2.,bbd/2.,bbd/2.);
    glScaled(0.1,0.1,0.1);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0,-100000);
    igl::opengl2::draw_beach_ball();
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPopMatrix();
  }

  // Mouse is always on top
  glDisable(GL_DEPTH_TEST);
  if(!panim.is_animating)
  {
    s.mouse.draw();
  }
  pop_object();
  pop_scene();

  igl::opengl::report_gl_error();

  TwDraw();
  glutSwapBuffers();
  if(canim.is_animating || panim.is_animating)
  {
    glutPostRedisplay();
  }
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
  glutPostRedisplay();
}

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  bool tw_using = TwEventMouseButtonGLUT(glutButton,glutState,mouse_x,mouse_y);
  const int mod = (glutButton <=2 ? glutGetModifiers() : 0);
  const bool option_down = mod & GLUT_ACTIVE_ALT;
  switch(glutButton)
  {
    case GLUT_RIGHT_BUTTON:
    case GLUT_LEFT_BUTTON:
    {
      push_scene();
      push_object();
      switch(glutState)
      {
        case 1:
        {
          // up
          const bool mouse_was_selecting = s.mouse.is_selecting();
          is_rotating = false;
          s.mouse.up(mouse_x,mouse_y);
          glutSetCursor(GLUT_CURSOR_INHERIT);
          if(mouse_was_selecting)
          {
            s.mouse.set_selection_from_last_drag(C,BE,P,RP);
            using namespace igl::opengl2;
            MouseController::VectorXb S;
            MouseController::propogate_to_descendants_if(
              s.mouse.selection(),P,S);
            MouseController::color_if(S,MAYA_SEA_GREEN,MAYA_VIOLET,s.colors);
          }
          break;
        }
        case 0:
          if(!tw_using)
          {
            down_x = mouse_x;
            down_y = mouse_y;
            if(option_down || glutButton==GLUT_RIGHT_BUTTON)
            {
              glutSetCursor(GLUT_CURSOR_CYCLE);
              // collect information for trackball
              is_rotating = true;
              down_camera = camera;
            }else
            {
              push_undo();
              s.mouse.down(mouse_x,mouse_y);
            }
          }
        break;
      }
      pop_object();
      pop_scene();
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
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  push_scene();
  push_object();
  if(is_rotating)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    Quaterniond q;
    switch(rotation_type)
    {
      case ROTATION_TYPE_IGL_TRACKBALL:
      {
        // Rotate according to trackball
        igl::trackball(
          width,
          height,
          2.0,
          down_camera.m_rotation_conj,
          down_x,
          down_y,
          mouse_x,
          mouse_y,
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
    camera.orbit(q.conjugate());
  }else if(s.mouse.drag(mouse_x,mouse_y))
  {
  }
  pop_object();
  pop_scene();
  glutPostRedisplay();
}

void init_relative()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  per_face_normals(V,F,N);
  const auto Vmax = V.colwise().maxCoeff();
  const auto Vmin = V.colwise().minCoeff();
  Vmid = 0.5*(Vmax + Vmin);
  bbd = (Vmax-Vmin).norm();
  camera.push_away(2);
}

void undo()
{
  using namespace std;
  if(!undo_stack.empty())
  {
    redo_stack.push(s);
    s = undo_stack.top();
    undo_stack.pop();
    s.mouse.reshape(width,height);
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
    s.mouse.reshape(width,height);
  }
}

bool save_pose()
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  string output_filename;
  next_filename(output_pose_prefix,4,".dmat",output_filename);
  MatrixXd T;
  forward_kinematics(C,BE,P,s.mouse.rotations(),T);
  if(writeDMAT(output_filename,T))
  {
    cout<<GREENGIN("Current pose written to "+output_filename+".")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Writing to "+output_filename+" failed.")<<endl;
    return false;
  }
}

bool save_weights()
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  if(writeDMAT(output_weights_filename,W))
  {
    cout<<GREENGIN("Current weights written to "+
      output_weights_filename+".")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Writing to "+output_weights_filename+" failed.")<<endl;
    return false;
  }
}

bool save_mesh()
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  MatrixXd T;
  forward_kinematics(C,BE,P,s.mouse.rotations(),T);
  MatrixXd U = M*T;
  string output_filename;
  next_filename(output_pose_prefix,4,".obj",output_filename);
  if(writeOBJ(output_filename,U,F))
  {
    cout<<GREENGIN("Current mesh written to "+
      output_filename+".")<<endl;
    return true;
  }else
  {
    cout<<REDRUM("Writing to "+output_filename+" failed.")<<endl;
    return false;
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
    case 'A':
    case 'a':
    {
      panim.is_animating = !panim.is_animating;
      panim.pose = s.mouse.rotations();
      panim.start_time = get_seconds();
      break;
    }
    case 'D':
    case 'd':
    {
      push_undo();
      s.mouse.clear_selection();
      break;
    }
    case 'R':
    {
      push_undo();
      s.mouse.reset_selected_rotations();
      break;
    }
    case 'r':
    {
      push_undo();
      s.mouse.reset_rotations();
      break;
    }
    case 'S':
    {
      save_mesh();
      break;
    }
    case 's':
    {
      save_pose();
      break;
    }
    case 'W':
    case 'w':
    {
      save_weights();
      break;
    }
    case 'z':
    case 'Z':
      is_rotating = false;
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
        snap_to_canonical_view_quat(camera.m_rotation_conj,1.0,q);
        camera.orbit(q.conjugate());
      }
      break;
    default:
      if(!TwEventKeyboardGLUT(key,mouse_x,mouse_y))
      {
        cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
      }
  }

  glutPostRedisplay();
}

bool clean(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  {
    MatrixXi _1;
    VectorXi _2,IM;
#ifdef VERBOSE
    cout<<"remesh_self_intersections"<<endl;
#endif
    igl::cgal::remesh_self_intersections(V,F,{},CV,CF,_1,_2,IM);
    for_each(CF.data(),CF.data()+CF.size(),[&IM](int & a){a=IM(a);});
    MatrixXd oldCV = CV;
    MatrixXi oldCF = CF;
    remove_unreferenced(oldCV,oldCF,CV,CF,IM);
  }
  MatrixXd TV;
  MatrixXi TT;
  {
    MatrixXi _1;
    // c  convex hull
    // Y  no boundary steiners
    // p  polygon input
#ifdef VERBOSE
    cout<<"tetrahedralize"<<endl;
#endif
    if(igl::tetgen::tetrahedralize(CV,CF,"cYpC",TV,TT,_1) != 0)
    {
      cout<<REDRUM("CDT failed.")<<endl;
      return false;
    }
  }
  {
    MatrixXd BC;
    barycenter(TV,TT,BC);
    VectorXd W;
#ifdef VERBOSE
    cout<<"winding_number"<<endl;
#endif
    winding_number(V,F,BC,W);
    W = W.array().abs();
    const double thresh = 0.5;
    const int count = (W.array()>thresh).cast<int>().sum();
    MatrixXi CT(count,TT.cols());
    int c = 0;
    for(int t = 0;t<TT.rows();t++)
    {
      if(W(t)>thresh)
      {
        CT.row(c++) = TT.row(t);
      }
    }
    assert(c==count);
    boundary_facets(CT,CF);
  }
  return true;
}

bool robust_weights(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & C,
  const Eigen::MatrixXi & BE,
  Eigen::MatrixXd & W)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  // clean mesh
  MatrixXd CV;
  MatrixXi CF;
  if(!clean(V,F,CV,CF))
  {
    return false;
  }
  MatrixXd TV;
  MatrixXi TT;
  // compute tet-mesh
  {
    MatrixXi _1;
#ifdef VERBOSE
    cout<<"mesh_with_skeleton"<<endl;
#endif
    if(!igl::tetgen::mesh_with_skeleton(CV,CF,C,{},BE,{},10,"pq1.5Y",TV,TT,_1))
    {
      cout<<REDRUM("tetgen failed.")<<endl;
      return false;
    }
  }
  // Finally, tetgen may have still included some insanely small tets.
  // Just ignore these during weight computation (and hope they don't isolate
  // any vertices).
  {
    const MatrixXi oldTT = TT;
    VectorXd vol;
    volume(TV,TT,vol);
    const double thresh = 1e-17;
    const int count = (vol.array()>thresh).cast<int>().sum();
    TT.resize(count,oldTT.cols());
    int c = 0;
    for(int t = 0;t<oldTT.rows();t++)
    {
      if(vol(t)>thresh)
      {
        TT.row(c++) = oldTT.row(t);
      }
    }
  }

  // compute weights
  VectorXi b;
  MatrixXd bc;
  if(!boundary_conditions(TV,TT,C,{},BE,{},b,bc))
  {
    cout<<REDRUM("boundary_conditions failed.")<<endl;
    return false;
  }
  // compute BBW
  // Default bbw data and flags
  igl::bbw::BBWData bbw_data;
  bbw_data.verbosity = 1;
#ifdef IGL_NO_MOSEK
  bbw_data.qp_solver = igl::bbw::QP_SOLVER_IGL_ACTIVE_SET;
  bbw_data.active_set_params.max_iter = 4;
#else
  bbw_data.mosek_data.douparam[MSK_DPAR_INTPNT_TOL_REL_GAP]=1e-14;
  bbw_data.qp_solver = igl::bbw::QP_SOLVER_MOSEK;
#endif
  // Weights matrix
  if(!igl::bbw::bbw(TV,TT,b,bc,bbw_data,W))
  {
    return false;
  }
  // Normalize weights to sum to one
  normalize_row_sums(W,W);
  // keep first #V weights
  W.conservativeResize(V.rows(),W.cols());
  return true;
}


int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  string filename = "../shared/cheburashka.off";
  string skel_filename = "../shared/cheburashka.tgf";
  string weights_filename = "";
  output_pose_prefix = "";
  switch(argc)
  {
    case 5:
      output_pose_prefix = argv[4];
      //fall through
    case 4:
      weights_filename = argv[3];
      //fall through
    case 3:
      skel_filename = argv[2];
      // Read and prepare mesh
      filename = argv[1];
      break;
    default:
      cerr<<"Usage:"<<endl<<"    ./example model.obj skeleton.tgf [weights.dmat] [pose-prefix]"<<endl;
      cout<<endl<<"Opening default rig..."<<endl;
  }

  // print key commands
  cout<<"[Click] and [drag]     Select a bone/Use onscreen widget to rotate bone."<<endl;
  cout<<"⌥ +[Click] and [drag]  Rotate secene."<<endl;
  cout<<"⌫                      Delete selected node(s) and incident bones."<<endl;
  cout<<"D,d                    Deselect all."<<endl;
  cout<<"R                      Reset selected rotation."<<endl;
  cout<<"r                      Reset all rotations."<<endl;
  cout<<"S                      Save current posed mesh."<<endl;
  cout<<"s                      Save current skeleton pose."<<endl;
  cout<<"W,w                    Save current weights."<<endl;
  cout<<"Z,z                    Snap to canonical view."<<endl;
  cout<<"⌘ Z                    Undo."<<endl;
  cout<<"⇧ ⌘ Z                  Redo."<<endl;
  cout<<"^C,ESC                 Exit (without saving)."<<endl;

  string dir,_1,_2,name;
  read_triangle_mesh(filename,V,F,dir,_1,_2,name);

  if(output_weights_filename.size() == 0)
  {
    output_weights_filename = dir+"/"+name+"-weights.dmat";
  }
  if(output_pose_prefix.size() == 0)
  {
    output_pose_prefix = dir+"/"+name+"-pose-";
  }

  {
    string output_filename;
    next_filename(output_pose_prefix,4,".dmat",output_filename);
    cout<<BLUEGIN("Output weights set to start with "<<
      output_weights_filename)<<endl;
    cout<<BLUEGIN("Output poses set to start with "<<output_filename)<<endl;
  }

  // Read in skeleton and precompute hierarchy
  readTGF(skel_filename,C,BE);
  // initialize mouse interface
  s.mouse.set_size(BE.rows());
  // Rigid parts (not used)
  colon<int>(0,BE.rows()-1,RP);
  assert(RP.size() == BE.rows());
  // Bone parents
  bone_parents(BE,P);
  if(weights_filename.size() == 0)
  {
    robust_weights(V,F,C,BE,W);
  }else
  {
    // Read in weights and precompute LBS matrix
    readDMAT(weights_filename,W);
  }
  lbs_matrix(V,W,M);

  init_relative();

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
    camera.m_rotation_conj.coeffs().data(), "open readonly=true");
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("RotationType",
    "igl_trackball,two-a...-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  rebar.TwAddVarRW("wireframe", TW_TYPE_BOOLCPP,&wireframe,"key=l");
  rebar.TwAddVarRW("centroid_is_visible", TW_TYPE_BOOLCPP,&centroid_is_visible,
    "keyIncr=C keyDecr=c label='centroid visible?'");
  TwType SkelStyleTypeTW = igl::anttweakbar::ReTwDefineEnumFromString("SkelStyleType",
    "3d,vector-graphics");
  rebar.TwAddVarRW("style",SkelStyleTypeTW,&skel_style,"");
  rebar.load(REBAR_NAME);

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT)/2.0);
  glutCreateWindow("skeleton-poser");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();

  return 0;
}

