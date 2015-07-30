#include <igl/Camera.h>
#include <igl/opengl/OpenGL_convenience.h>
#include <igl/PI.h>
#include <igl/STR.h>
#include <igl/arap.h>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/get_seconds.h>
#include <igl/harmonic.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/list_to_matrix.h>
#include <igl/massmatrix.h>
#include <igl/material_colors.h>
#include <igl/matlab_format.h>
#include <igl/normalize_row_lengths.h>
#include <igl/partition.h>
#include <igl/pathinfo.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/quat_to_mat.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/readWRL.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/snap_to_fixed_up.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/writeDMAT.h>
#include <igl/anttweakbar/ReAntTweakBar.h>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include <Eigen/Core>

#include <vector>
#include <iostream>
#include <algorithm>

struct State
{
  igl::Camera camera;
} s;
enum RotationType
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type;
bool is_rotating = false;
int down_x,down_y;
igl::Camera down_camera;
bool is_animating = false;
double animation_start_time = 0;
double ANIMATION_DURATION = 0.5;
Eigen::Quaterniond animation_from_quat;
Eigen::Quaterniond animation_to_quat;
// Use vector for range-based `for`
std::vector<State> undo_stack;
std::vector<State> redo_stack;
bool paused = false;
double t = 0;
double pause_time = 0.0;

void push_undo()
{
  undo_stack.push_back(s);
  // Clear
  redo_stack = std::vector<State>();
}

void undo()
{
  using namespace std;
  if(!undo_stack.empty())
  {
    redo_stack.push_back(s);
    s = undo_stack.front();
    undo_stack.pop_back();
  }
}

void redo()
{
  using namespace std;
  if(!redo_stack.empty())
  {
    undo_stack.push_back(s);
    s = redo_stack.front();
    redo_stack.pop_back();
  }
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

// Width and height of window
int width,height;
// Position of light
float light_pos[4] = {0.1,0.1,-0.9,0};
// Vertex positions, normals, colors and centroid
Eigen::MatrixXd V,U,N,C,mid;
Eigen::VectorXi S;
igl::ARAPData arap_data;
Eigen::MatrixXi F,T;
int selected_col = 0;
// Faces
// Bounding box diagonal length
double bbd;
int tot_num_samples = 0;
#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar; // Pointer to the tweak bar

int num_in_selection(const Eigen::VectorXi & S)
{
  int count = 0;
  for(int v = 0;v<S.rows(); v++)
  {
    if(S(v) >= 0)
    {
      count++;
    }
  }
  return count;
}

bool init_arap()
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  VectorXi b(num_in_selection(S));
  assert(S.rows() == V.rows());
  C.resize(S.rows(),3);
  MatrixXd bc = MatrixXd::Zero(b.size(),S.maxCoeff()+1);
  MatrixXi * Ele;
  if(T.rows()>0)
  {
    Ele = &T;
  }else
  {
    Ele = &F;
  }
  // get b from S
  {
    int bi = 0;
    for(int v = 0;v<S.rows(); v++)
    {
      if(S(v) >= 0)
      {
        b(bi) = v;
        bc(bi,S(v)) = 1;
        bi++;
        switch(S(v))
        {
          case 0:
            C.row(v) = RowVector3d(0.039,0.31,1);
            break;
          case 1:
            C.row(v) = RowVector3d(1,0.41,0.70);
            break;
          default:
            C.row(v) = RowVector3d(0.4,0.8,0.3);
            break;
        }
      }else
      {
        C.row(v) = RowVector3d(
          GOLD_DIFFUSE[0],
          GOLD_DIFFUSE[1],
          GOLD_DIFFUSE[2]);
      }
    }
  }
  // Store current mesh
  U = V;
  VectorXi _S;
  VectorXd _D;
  MatrixXd W;
  if(!harmonic(V,*Ele,b,bc,1,W))
  {
    return false;
  }
  arap_data.with_dynamics = true;
  arap_data.h = 0.05;
  //arap_data.max_iter = 100;
  //partition(W,100,arap_data.G,_S,_D);
  return arap_precomputation(V,*Ele,V.cols(),b,arap_data);
}

bool update_arap()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  MatrixXd bc(num_in_selection(S),V.cols());
  // get b from S
  {
    if(!paused)
    {
      t = get_seconds()-pause_time;
    }
    int bi = 0;
    for(int v = 0;v<S.rows(); v++)
    {
      if(S(v) >= 0)
      {
        bc.row(bi) = V.row(v);
        switch(S(v))
        {
          case 0:
          {
            const double r = mid(0)*0.25;
            bc(bi,0) += r*cos(0.5*t*2.*PI);
            bc(bi,1) -= r+r*sin(0.5*t*2.*PI);
            break;
          }
          case 1:
          {
            const double r = mid(1)*0.15;
            bc(bi,1) += r+r*cos(0.15*t*2.*PI);
            bc(bi,2) -= r*sin(0.15*t*2.*PI);

            //// Pull-up
            //bc(bi,0) += 0.42;//mid(0)*0.5;
            //bc(bi,1) += 0.55;//mid(0)*0.5;
            //// Bend
            //Vector3d t(-1,0,0);
            //Quaterniond q(AngleAxisd(PI/1.5,Vector3d(0,1.0,0.1).normalized()));
            //const Vector3d a = bc.row(bi);
            //bc.row(bi) = (q*(a-t) + t) + Vector3d(1.5,0.1,0.9);


            break;
          }
          default:
          break;
        }
        bi++;
      }
    }
  }
  if(!arap_solve(bc,arap_data,U))
  {
    cerr<<"arap_solve failed."<<endl;
    return false;
  }
  per_face_normals(U,F,N);
  return true;
}

void reshape(int width,int height)
{
  using namespace std;
  // Save width and height
  ::width = width;
  ::height = height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,width,height);
  // Send the new window size to AntTweakBar
  TwWindowSize(width, height);
  // Set aspect for all cameras
  s.camera.m_aspect = (double)width/(double)height;
  for(auto & s : undo_stack)
  {
    s.camera.m_aspect = (double)width/(double)height;
  }
  for(auto & s : redo_stack)
  {
    s.camera.m_aspect = (double)width/(double)height;
  }
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

void pop_scene()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void pop_object()
{
  glPopMatrix();
}

// Scale and shift for object
void push_object()
{
  glPushMatrix();
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-mid(0,0),-mid(0,1),-mid(0,2));
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  //glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  float amb[4];
  amb[0] = amb[1] = amb[2] = 0;
  amb[3] = 1.0;
  float diff[4] = {0.0,0.0,0.0,0.0};
  diff[0] = diff[1] = diff[2] = (1.0 - 0/0.4);;
  diff[3] = 1.0;
  float zeros[4] = {0.0,0.0,0.0,0.0};
  float pos[4];
  copy(light_pos,light_pos+4,pos);
  glLightfv(GL_LIGHT0,GL_AMBIENT,amb);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,diff);
  glLightfv(GL_LIGHT0,GL_SPECULAR,zeros);
  glLightfv(GL_LIGHT0,GL_POSITION,pos);
  pos[0] *= -1;
  pos[1] *= -1;
  pos[2] *= -1;
  glLightfv(GL_LIGHT1,GL_AMBIENT,amb);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,diff);
  glLightfv(GL_LIGHT1,GL_SPECULAR,zeros);
  glLightfv(GL_LIGHT1,GL_POSITION,pos);
}

//const float back[4] = {30.0/255.0,30.0/255.0,50.0/255.0,0};
const float back[4] = {255.0/255.0,255.0/255.0,255.0/255.0,0};
void display()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  // Update
  update_arap();

  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(is_animating)
  {
    double t = (get_seconds() - animation_start_time)/ANIMATION_DURATION;
    if(t > 1)
    {
      t = 1;
      is_animating = false;
    }
    const Quaterniond q = animation_from_quat.slerp(t,animation_to_quat).normalized();
    s.camera.orbit(q.conjugate());
  }

  glDisable(GL_LIGHTING);
  lights();
  push_scene();
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_NORMALIZE);
  push_object();

  // Draw the model
  // Set material properties
  //glDisable(GL_COLOR_MATERIAL);
  //glMaterialfv(GL_FRONT, GL_AMBIENT,  GOLD_AMBIENT);
  //glMaterialfv(GL_FRONT, GL_DIFFUSE,  GOLD_DIFFUSE  );
  //glMaterialfv(GL_FRONT, GL_SPECULAR, GOLD_SPECULAR);
  //glMaterialf (GL_FRONT, GL_SHININESS, 128);
  //glMaterialfv(GL_BACK, GL_AMBIENT,  SILVER_AMBIENT);
  //glMaterialfv(GL_BACK, GL_DIFFUSE,  FAST_GREEN_DIFFUSE  );
  //glMaterialfv(GL_BACK, GL_SPECULAR, SILVER_SPECULAR);
  //glMaterialf (GL_BACK, GL_SHININESS, 128);
  glEnable(GL_COLOR_MATERIAL);

  igl::opengl2::draw_mesh(U,F,N,C);
  glDisable(GL_COLOR_MATERIAL);

  pop_object();


  // Draw a nice floor
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-mid(1));
  glTranslated(0,floor_offset,0);
  //const float GREY[4] = {0.5,0.5,0.6,1.0};
  //const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};

  //draw_floor(GREY,DARK_GREY);
  igl::opengl2::draw_floor();
  glPopMatrix();

  pop_scene();

  igl::opengl::report_gl_error();

  TwDraw();
  glutSwapBuffers();
  //if(is_animating)
  //{
    glutPostRedisplay();
  //}
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
          glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
          is_rotating = false;
          break;
        case 0:
          // down
          if(!tw_using)
          {
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
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
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
  }else
  {
    TwEventMouseMotionGLUT(mouse_x, mouse_y);
  }
  glutPostRedisplay();
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace igl;
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(0);
    case ' ':
      {
        static double pause_start,pause_stop;
        paused = !paused;
        if(paused)
        {
          pause_start = get_seconds();
        }else
        {
          pause_stop = get_seconds();
          pause_time += (pause_stop-pause_start);
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
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  // init mesh
  string filename = "../shared/decimated-knight.obj";
  string sfilename = "../shared/decimated-knight-selection.dmat";
  //string filename = "../shared/decimated-knight.mesh";
  //string sfilename = "../shared/decimated-knight-1-selection.dmat";
  if(argc < 3)
  {
    cerr<<"Usage:"<<endl<<"    ./example input.obj selection.dmat"<<endl;
    cout<<endl<<"Opening default mesh..."<<endl;
  }else
  {
    // Read and prepare mesh
    filename = argv[1];
    sfilename = argv[2];
  }
  string d,b,ext,f;
  pathinfo(filename,d,b,ext,f);
  // Convert extension to lower case
  transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  vector<vector<double > > vV,vN,vTC;
  vector<vector<int > > vF,vTF,vFN;
  // Convert extension to lower case
  if(ext == "obj")
  {
    if(!igl::readOBJ(filename,vV,vTC,vN,vF,vTF,vFN))
    {
      return 1;
    }
  }else if(ext == "mesh")
  {
    if(!igl::readMESH(filename,V,T,F))
    {
      return 1;
    }
  }else
  {
    return 1;
  }
  if(vV.size() > 0)
  {
    if(!list_to_matrix(vV,V))
    {
      cerr<<"Bad V"<<endl;
      return 1;
    }
    polygon_mesh_to_triangle_mesh(vF,F);
  }
  per_face_normals(V,F,N);

  if(!readDMAT(sfilename,S))
  {
    return 1;
  }

  // Compute normals, centroid, colors, bounding box diagonal
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  bbd = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();

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
  s.camera.push_away(3);
  s.camera.dolly_zoom(25-s.camera.m_angle);
  TwType RotationTypeTW = igl::anttweakbar::ReTwDefineEnumFromString(
     "RotationType",
    "igl_trackball,two-a...-fixed-up");
  rebar.TwAddVarCB( "rotation_type", RotationTypeTW,
    set_rotation_type,get_rotation_type,NULL,"keyIncr=] keyDecr=[");
  rebar.load(REBAR_NAME);

  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("colored-mesh");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc(
    [](int x, int y)
    {
      TwEventMouseMotionGLUT(x,y);
      glutPostRedisplay();
    });
  static std::function<void(int)> timer_bounce;
  auto timer = [] (int ms) {
    timer_bounce(ms);
  };
  timer_bounce = [&] (int ms) {
    glutTimerFunc(ms, timer, ms);
    glutPostRedisplay();
  };
  glutTimerFunc(500, timer, 500);

  if(!init_arap())
  {
    cerr<<"Initializing arap failed."<<endl;
    return 1;
  }

  glutMainLoop();
  return 0;
}
