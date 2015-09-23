#include <igl/opengl/OpenGL_convenience.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/normalize_row_lengths.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/quat_to_mat.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/readMESH.h>
#include <igl/readWRL.h>
#include <igl/trackball.h>
#include <igl/list_to_matrix.h>
#include <igl/polygon_mesh_to_triangle_mesh.h>
#include <igl/material_colors.h>
#include <igl/barycenter.h>
#include <igl/matlab_format.h>
#include <igl/anttweakbar/ReAntTweakBar.h>
#include <igl/pathinfo.h>
#include <igl/embree/EmbreeIntersector.h>
#include <igl/embree/ambient_occlusion.h>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include <Eigen/Core>

#include <vector>
#include <iostream>
#include <algorithm>

// Width and height of window
int width,height;
// Rotation of scene
float scene_rot[4] = {0,0,0,1};
// information at mouse down
float down_scene_rot[4] = {0,0,0,1};
bool trackball_on = false;
int down_mouse_x,down_mouse_y;
// Position of light
float light_pos[4] = {0.1,0.1,-0.9,0};
// Vertex positions, normals, colors and centroid
Eigen::MatrixXd V,N,C,mid;
// Faces
Eigen::MatrixXi F;
// Bounding box diagonal length
double bbd;
igl::embree::EmbreeIntersector ei;
// Running ambient occlusion
Eigen::VectorXd S;
int tot_num_samples = 0;
#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar; // Pointer to the tweak bar
bool lights_on = true;
Eigen::Vector4f color(0.4,0.8,0.3,1.0);
double ao_factor = 1.0;
bool ao_normalize = false;
bool ao_on = true;
double light_intensity = 1.0;

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
}

// Set up igl::opengl2::projection and model view of scene
void push_scene()
{
  using namespace igl;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45,(double)width/(double)height,1e-2,100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0,0,3,0,0,0,0,1,0);
  glPushMatrix();
  float mat[4*4];
  quat_to_mat(scene_rot,mat);
  glMultMatrixf(mat);
}

void pop_scene()
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

void pop_object()
{
  glPopMatrix();
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  float amb[4];
  amb[0] = amb[1] = amb[2] = light_intensity;
  amb[3] = 1.0;
  float diff[4] = {0.0,0.0,0.0,0.0};
  diff[0] = diff[1] = diff[2] = (1.0 - light_intensity/0.4);;
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

const float back[4] = {30.0/255.0,30.0/255.0,50.0/255.0,0};
void display()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  if(!trackball_on && tot_num_samples < 10000)
  {
    if(S.size() == 0)
    {
      S.resize(V.rows());
      S.setZero();
    }
    VectorXd Si;
    const int num_samples = 20;
    igl::embree::ambient_occlusion(ei,V,N,num_samples,Si);
    S *= (double)tot_num_samples;
    S += Si*(double)num_samples;
    tot_num_samples += num_samples;
    S /= (double)tot_num_samples;
  }

  // Convert to 1-intensity
  C.conservativeResize(S.rows(),3);
  if(ao_on)
  {
    C<<S,S,S;
    C.array() = (1.0-ao_factor*C.array());
  }else
  {
    C.setConstant(1.0);
  }
  if(ao_normalize)
  {
    C.col(0) *= ((double)C.rows())/C.col(0).sum();
    C.col(1) *= ((double)C.rows())/C.col(1).sum();
    C.col(2) *= ((double)C.rows())/C.col(2).sum();
  }
  C.col(0) *= color(0);
  C.col(1) *= color(1);
  C.col(2) *= color(2);

  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // All smooth points
  glEnable( GL_POINT_SMOOTH );

  glDisable(GL_LIGHTING);
  if(lights_on)
  {
    lights();
  }
  push_scene();
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  push_object();

  // Draw the model
  // Set material properties
  glEnable(GL_COLOR_MATERIAL);
  igl::opengl2::draw_mesh(V,F,N,C);

  pop_object();

  // Draw a nice floor
  glPushMatrix();
  const double floor_offset =
    -2./bbd*(V.col(1).maxCoeff()-mid(1));
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

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  bool tw_using = TwEventMouseButtonGLUT(glutButton,glutState,mouse_x,mouse_y);
  switch(glutState)
  {
    case 1:
      // up
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
      trackball_on = false;
      break;
    case 0:
      // down
      if(!tw_using)
      {
        glutSetCursor(GLUT_CURSOR_CYCLE);
        // collect information for trackball
        trackball_on = true;
        copy(scene_rot,scene_rot+4,down_scene_rot);
        down_mouse_x = mouse_x;
        down_mouse_y = mouse_y;
      }
    break;
  }
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;


  if(trackball_on)
  {
    // Rotate according to trackball
    trackball<float>(
      width,
      height,
      2,
      down_scene_rot,
      down_mouse_x,
      down_mouse_y,
      mouse_x,
      mouse_y,
      scene_rot);
  }else
  {
    TwEventMouseMotionGLUT(mouse_x, mouse_y);
  }
}


void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  switch(key)
  {
    // ESC
    case char(27):
      rebar.save(REBAR_NAME);
    // ^C
    case char(3):
      exit(0);
    default:
      if(!TwEventKeyboardGLUT(key,mouse_x,mouse_y))
      {
        cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
      }
  }

}


int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  // init mesh
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

  // Compute normals, centroid, colors, bounding box diagonal
  per_vertex_normals(V,F,N);
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  bbd = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();

  // Init embree
  ei.init(V.cast<float>(),F.cast<int>());

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
  rebar.TwAddVarRW("scene_rot", TW_TYPE_QUAT4F, &scene_rot, "");
  rebar.TwAddVarRW("lights_on", TW_TYPE_BOOLCPP, &lights_on, "key=l");
  rebar.TwAddVarRW("color", TW_TYPE_COLOR4F, color.data(), "colormode=hls");
  rebar.TwAddVarRW("ao_factor", TW_TYPE_DOUBLE, &ao_factor, "min=0 max=1 step=0.2 keyIncr=] keyDecr=[ ");
  rebar.TwAddVarRW("ao_normalize", TW_TYPE_BOOLCPP, &ao_normalize, "key=n");
  rebar.TwAddVarRW("ao_on", TW_TYPE_BOOLCPP, &ao_on, "key=a");
  rebar.TwAddVarRW("light_intensity", TW_TYPE_DOUBLE, &light_intensity, "min=0 max=0.4 step=0.1 keyIncr=} keyDecr={ ");
  rebar.load(REBAR_NAME);

  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("ambient-occlusion");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
  glutMainLoop();
  return 0;
}
