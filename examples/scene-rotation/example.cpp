// Small GLUT application to test different scene rotation paradigms 
//

#include "trackball.h"

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

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <GLUT/glut.h>

#include <Carbon/Carbon.h>

#include <string>
#include <vector>
#include <stack>
#include <iostream>


Eigen::MatrixXd V,N;
Eigen::VectorXd Vmid,Vmin,Vmax;
double bbd = 1.0;
Eigen::MatrixXi F;
struct State
{
  igl::Camera camera;
} s;

// See README for descriptions
enum ROTATION_TYPE
{
  ROTATION_TYPE_IGL_TRACKBALL = 0,
  ROTATION_TYPE_BELL_TRACKBALL = 1,
  ROTATION_TYPE_TWO_AXIS_VALUATOR = 2,
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 3,
  NUM_ROTATION_TYPES = 4,
} rotation_type;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool is_rotating = false;
int down_x,down_y;
igl::Camera down_camera;

int width,height;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

#define REBAR_NAME "temp.rbr"
igl::ReTwBar rebar;

// No-op setter, does nothing
void TW_CALL no_op(const void * /*value*/, void * /*clientData*/)
{
}

void TW_CALL get_camera_rotation(void * value, void *clientData)
{
  using namespace std;
  // case current value to double
  double * quat = (double *)(value);
  std::copy(s.camera.rotation,s.camera.rotation+4,quat);
}

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

  
  draw_mesh(V,F,N);
  pop_object();

  // Draw a nice floor
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-Vmid(1));
  glTranslated(0,floor_offset,0);
  const float GREY[4] = {0.5,0.5,0.6,1.0};
  const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
  draw_floor(GREY,DARK_GREY);
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
  push_undo();
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
      case ROTATION_TYPE_BELL_TRACKBALL:
      {
        float down_quaternion[4];
        copy(down_camera.rotation,down_camera.rotation+4,down_quaternion);
        float new_quaternion[4];
        
        const float center_x = ((float)width)/2.0;
        const float center_y = ((float)height)/2.0;
        const double speed = 2.0f;
        const float half_width =   ((float)width)/speed;
        const float half_height = ((float)height)/speed;
        
        ::trackball(new_quaternion,
          (float)(center_x-down_x)/half_width,
          (float)(down_y-center_y)/half_height,
          (float)(center_x-mouse_x)/half_width,
          (float)(mouse_y-center_y)/half_height);
        // I think we need to do this because we have z pointing out of the
        // screen rather than into the screen
        new_quaternion[2] = -new_quaternion[2];
        float float_quat[4];
        add_quats(down_quaternion,new_quaternion,float_quat);
        copy(float_quat,float_quat+4,s.camera.rotation);
        break;
      }
      case ROTATION_TYPE_TWO_AXIS_VALUATOR:
      {
        Quaterniond down_q;
        copy(down_camera.rotation,down_camera.rotation+4,down_q.coeffs().data());
        Vector3d axis(mouse_y-down_y,mouse_x-down_x,0);
        const double speed = 2.0;
        if(axis.norm() != 0)
        {
          Quaterniond q;
          q = 
            Quaterniond(
              AngleAxisd(
                M_PI*axis.norm()/(double)width*speed/2.0,
                axis.normalized())) * down_q;
          q.normalize();
          copy(q.coeffs().data(),q.coeffs().data()+4,s.camera.rotation);
        }
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


KeyMap keyStates ;
bool IS_KEYDOWN( uint16_t vKey )
{
  uint8_t index = vKey / 32 ;
  uint8_t shift = vKey % 32 ;
  return keyStates[index].bigEndianValue & (1 << shift) ;
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
  GetKeys(keyStates);
  const bool command_down = IS_KEYDOWN(kVK_Command);
  const bool shift_down = IS_KEYDOWN(kVK_Shift);
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
  string filename = "../shared/cheburashka.obj";
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
  //    boundary_faces(T,F);
  //  }
  }
  if(vV.size() > 0)
  {
    if(!list_to_matrix(vV,V))
    {
      return 1;
    }
    triangulate(vF,F);
  }

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
  rebar.TwAddVarCB("camera_rotation", TW_TYPE_QUAT4D, no_op,get_camera_rotation, NULL, "open readonly=true");
  TwEnumVal RotationTypesEV[NUM_ROTATION_TYPES] = 
  {
    {ROTATION_TYPE_IGL_TRACKBALL,"igl trackball"},
    {ROTATION_TYPE_BELL_TRACKBALL,"bell trackball"},
    {ROTATION_TYPE_TWO_AXIS_VALUATOR,"two axis valuator"},
    {ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP,"two a... fixed up"},
  };
  TwType RotationTypeTW = 
    ReTwDefineEnum(
        "RotationType", 
        RotationTypesEV, 
        NUM_ROTATION_TYPES);
  rebar.TwAddVarRW( "rotation_type", RotationTypeTW, &rotation_type,"keyIncr=] keyDecr=[");
  rebar.load(REBAR_NAME);

  // Init antweakbar
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
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
