// Small GLUT application to test different scene rotation paradigms 
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
#include <igl/PI.h>
#include <igl/render_to_tga.h>
#include <igl/STR.h>
#define IGL_HEADER_ONLY
#include <igl/lens_flare.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <GLUT/glut.h>

#include <Carbon/Carbon.h>

#include <string>
#include <vector>
#include <iomanip>
#include <stack>
#include <iostream>

bool eyes_visible = true;
double x=6,y=232,z=61;

std::vector<igl::Flare> flares;
std::vector<GLuint> shine_ids;
std::vector<GLuint> flare_ids;
int shine_tic = 0;

GLuint list_id = 0;
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
  ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
  NUM_ROTATION_TYPES = 2,
} rotation_type;

std::stack<State> undo_stack;
std::stack<State> redo_stack;

bool is_rotating = false;
int down_x,down_y;
igl::Camera down_camera;
bool render_to_tga_on_next = false;
int render_count = 0;

int width,height;
Eigen::Vector4f light_pos(-0.1,-0.1,0.9,0);

#define REBAR_NAME "temp.rbr"
igl::ReTwBar rebar;

void TW_CALL set_camera_rotation(const void * value, void *clientData)
{
  using namespace std;
  // case current value to double
  const double * quat = (const double *)(value);
  std::copy(quat,quat+4,s.camera.rotation);
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
  float GREY[4] =  {0.2,0.2,0.2,1.};
  float BLACK[4] =  {0.,0.,0.,1.};
  Vector4f pos = light_pos;
  glLightfv(GL_LIGHT0,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT0,GL_SPECULAR,GREY);
  glLightfv(GL_LIGHT0,GL_POSITION,pos.data());
  pos(0) *= -1;
  pos(1) *= -1;
  pos(2) *= -1;
  glLightfv(GL_LIGHT1,GL_AMBIENT,GREY);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,WHITE);
  glLightfv(GL_LIGHT1,GL_SPECULAR,BLACK);
  glLightfv(GL_LIGHT1,GL_POSITION,pos.data());
}

void init_flare()
{
}

void draw_flare()
{
  using namespace igl;
  using namespace Eigen;
  glPushMatrix();
  glScaled(bbd*0.5,bbd*0.5,bbd*0.5);
  glScaled(0.2,0.2,0.2);
  Vector3f light(0,0,0);
  lens_flare_draw(flares,shine_ids,flare_ids,light,1.0,shine_tic);
  glPopMatrix();
}

void draw_eyes()
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
#define NUM_LEDS 2
  Vector3d LED_pos[NUM_LEDS];
  LED_pos[0] = Vector3d( x,y,z);
  LED_pos[1] = Vector3d(-x,y,z);
  enum LEDMethod
  {
    LED_METHOD_COLORED_CIRCLE = 0,
    LED_METHOD_OUTLINED_CIRCLE = 1,
    LED_METHOD_TEXTURE_FLARE = 2
  } method = LED_METHOD_COLORED_CIRCLE;

  
  for(int l = 0;l<NUM_LEDS;l++)
  {
    glPushMatrix();
    glTranslated(LED_pos[l](0), LED_pos[l](1), LED_pos[l](2));
    const double r = 2.0;
    const float color[4] = {1,0,0,1};
    glScaled(r,r,r);
    switch(method)
    {
      case LED_METHOD_COLORED_CIRCLE:
      {
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT,GL_AMBIENT);
        glColor4fv(color);
        glBegin(GL_TRIANGLE_FAN);
        glVertex3d(0,0,0);
        for(double theta = 0;theta<2.*PI;theta+=2.*PI/15.)
        {
          glVertex3d(cos(theta),sin(theta),0);
        }
        glEnd();
        break;
      }
      case LED_METHOD_OUTLINED_CIRCLE:
      {
        glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_LIGHTING);
        glColorMaterial(GL_FRONT,GL_DIFFUSE);
        glBegin(GL_TRIANGLE_FAN);
        glColor4fv(color);
        glVertex3d(0,0,0);
        glColor4f(0,0,0,1);
        for(double theta = 0;theta<2.*PI;theta+=2.*PI/15.)
        {
          glVertex3d(cos(theta),sin(theta),0);
        }
        glEnd();
        break;
      }
      case LED_METHOD_TEXTURE_FLARE:
      {
        draw_flare();
        break;
      }
    }
    glPopMatrix();
  }
}

void display()
{
  using namespace igl;
  using namespace std;
  glClearColor(0.03,0.03,0.04,0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_NORMALIZE);
  lights();
  push_scene();


  
  if(list_id == 0)
  {
    list_id = glGenLists(1);
    glNewList(list_id,GL_COMPILE);

    push_object();
    // Set material properties
    glDisable(GL_COLOR_MATERIAL);
    const float NEAR_BLACK[4] = {0.1,0.1,0.1,1.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  NEAR_BLACK);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  MIDNIGHT_BLUE_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SILVER_SPECULAR);
    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 128);
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

    glEndList();
  }
  glCallList(list_id);


  push_object();
  if(eyes_visible)
  {
    draw_eyes();
  }
  pop_object();

  pop_scene();

  report_gl_error();

  if(render_to_tga_on_next)
  {
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    render_to_tga(
      STR("./"<< "flare-eyes-" << setw(4) << setfill('0') << render_count++ << ".tga"),
      viewport[2],viewport[3],true);
    //render_to_tga_on_next = false;
  }

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
  using namespace std;
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
    case ' ':
      render_to_tga_on_next = !render_to_tga_on_next;
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
  string filename = "../shared/beast.obj";
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
  rebar.TwAddVarCB("camera_rotation", TW_TYPE_QUAT4D, set_camera_rotation,get_camera_rotation, NULL, "open");
  TwEnumVal RotationTypesEV[NUM_ROTATION_TYPES] = 
  {
    {ROTATION_TYPE_IGL_TRACKBALL,"igl trackball"},
    {ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP,"two a... fixed up"},
  };
  TwType RotationTypeTW = 
    ReTwDefineEnum(
        "RotationType", 
        RotationTypesEV, 
        NUM_ROTATION_TYPES);
  rebar.TwAddVarRW( "rotation_type", RotationTypeTW, &rotation_type,"keyIncr=] keyDecr=[");
  rebar.TwAddVarRW( "x",TW_TYPE_DOUBLE, &x,"");
  rebar.TwAddVarRW( "y",TW_TYPE_DOUBLE, &y,"");
  rebar.TwAddVarRW( "z",TW_TYPE_DOUBLE, &z,"");
  rebar.TwAddVarRW( "eyes_visible",TW_TYPE_BOOLCPP, &eyes_visible,"key=e");
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
  // Init flares
  lens_flare_load_textures(shine_ids,flare_ids);
  const float RED[3] = {1,0,0};
  const float GREEN[3] = {0,1,0};
  const float BLUE[3] = {0,0,1};
  //lens_flare_create(RED,GREEN,BLUE,flares);
  flares.resize(4);
  flares[0] = Flare(-1, 1.0f, 1.*0.1f,  RED, 1.0);
  flares[1] = Flare(-1, 1.0f, 1.*0.15f, GREEN, 1.0);
  flares[2] = Flare(-1, 1.0f, 1.*0.35f, BLUE, 1.0);
  flares[3] = Flare( 2, 1.0f, 1.*0.1f, BLUE, 0.4);

  glutMainLoop();

  return 0;
}
