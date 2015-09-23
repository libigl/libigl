#include <igl/EPS.h>
#include <igl/opengl/OpenGL_convenience.h>
#include <igl/colon.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/jet.h>
#include <igl/material_colors.h>
#include <igl/matlab_format.h>
#include <igl/normalize_row_lengths.h>
#include <igl/per_face_normals.h>
#include <igl/quat_to_mat.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/slice.h>
#include <igl/opengl2/sort_triangles.h>
#include <igl/trackball.h>
#include <igl/opengl2/unproject.h>
#include <igl/anttweakbar/ReAntTweakBar.h>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <Eigen/Core>

#include <vector>
#include <iostream>

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
Eigen::MatrixXd V,N,sorted_N,C,mid;
Eigen::VectorXi I;
// Bounding box diagonal length
double bbd;
// Faces
Eigen::MatrixXi F,sorted_F;
// alpha
double alpha = 0.2;

#define REBAR_NAME "temp.rbr"
igl::anttweakbar::ReTwBar rebar; // Pointer to the tweak bar

void reshape(int width,int height)
{
  using namespace std;
  // Save width and height
  ::width = width;
  ::height = height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,width,height);
}

// Set up double-sided lights
void lights()
{
  using namespace std;
  glEnable(GL_LIGHTING);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  float ones[4] = {1.0,1.0,1.0,1.0};
  float zeros[4] = {0.0,0.0,0.0,0.0};
  float pos[4];
  copy(light_pos,light_pos+4,pos);
  glLightfv(GL_LIGHT0,GL_AMBIENT,zeros);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,ones);
  glLightfv(GL_LIGHT0,GL_SPECULAR,zeros);
  glLightfv(GL_LIGHT0,GL_POSITION,pos);
  pos[0] *= -1;
  pos[1] *= -1;
  pos[2] *= -1;
  glLightfv(GL_LIGHT1,GL_AMBIENT,zeros);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,ones);
  glLightfv(GL_LIGHT1,GL_SPECULAR,zeros);
  glLightfv(GL_LIGHT1,GL_POSITION,pos);
}

// Set up igl::opengl2::projection and model view of scene
void push_scene()
{
  using namespace igl;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(45,(double)width/(double)height,1e-2,20);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(0,0,3,0,0,0,0,1,0);
  glPushMatrix();
  float mat[4*4];
  quat_to_mat(scene_rot,mat);
  glMultMatrixf(mat);
}

void pop_scene()
{
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
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

void init_C(const Eigen::MatrixXi & F)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  C.resize(F.rows(),3);
  for(int c = 0;c<F.rows();c++)
  {
    C(c,0) = C(c,1) = C(c,2) = (double)c/F.rows();
    //jet((double)c/F.rows(),C(c,0),C(c,1),C(c,2));
  }
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
      glutSetCursor(GLUT_CURSOR_INHERIT);
      trackball_on = false;
      // sort!
      push_scene();
      push_object();
      igl::opengl2::sort_triangles(V,F,sorted_F,I);
      slice(N,I,1,sorted_N);
      init_C(I);
      pop_object();
      pop_scene();
      break;
    case 0:
      // down
      glutSetCursor(GLUT_CURSOR_CYCLE);
      // collect information for trackball
      trackball_on = true;
      copy(scene_rot,scene_rot+4,down_scene_rot);
      down_mouse_x = mouse_x;
      down_mouse_y = mouse_y;
    break;
  }
}

const float BG[4] = {190.0/255.0,190.0/255.0,190.0/255.0,0};
void display()
{
  if(sorted_F.size() != F.size())
  {
    mouse(0,1,-1,-1);
  }
  //using namespace Eigen;
  using namespace igl;
  using namespace std;
  glClearColor(BG[0],BG[1],BG[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // All smooth points
  glEnable( GL_POINT_SMOOTH );

  lights();
  push_scene();
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glEnable(GL_CULL_FACE);
  //glCullFace(GL_BACK);

  // Draw a nice floor
  push_object();
  glEnable(GL_LIGHTING);
  glTranslated(0,V.col(1).minCoeff(),0);
  glScaled(2*bbd,2*bbd,2*bbd);
  glTranslated(0,-1000*FLOAT_EPS,0);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  igl::opengl2::draw_floor();
  pop_object();

  glDisable(GL_CULL_FACE);
  if(trackball_on)
  {
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }else
  {
    float front[4];
    copy(GOLD_DIFFUSE,GOLD_DIFFUSE+3,front);
    float back[4];
    //copy(FAST_GREEN_DIFFUSE,FAST_GREEN_DIFFUSE+3,back);
    copy(GOLD_DIFFUSE,GOLD_DIFFUSE+3,back);
    float amb[4];
    copy(SILVER_AMBIENT,SILVER_AMBIENT+3,amb);
    float spec[4];
    copy(SILVER_SPECULAR,SILVER_SPECULAR+3,spec);
    front[3] = back[3] = amb[3] = spec[3] = alpha;
    // Set material properties
    glDisable(GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT, GL_AMBIENT,  amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  front);
    glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
    glMaterialf (GL_FRONT, GL_SHININESS, 128);
    glMaterialfv(GL_BACK, GL_AMBIENT,  amb);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  back);
    glMaterialfv(GL_BACK, GL_SPECULAR, spec);
    glMaterialf (GL_BACK, GL_SHININESS, 128);
    glEnable(GL_LIGHTING);
  }

  push_object();
  // Draw the model
  igl::opengl2::draw_mesh(V,sorted_F,sorted_N,C);
  pop_object();

  pop_scene();

  igl::opengl::report_gl_error();
  glutSwapBuffers();
  glutPostRedisplay();
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
  }
}

void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  switch(key)
  {
    // Ctrl-c and esc exit
    case char(3):
    case char(27):
      rebar.save(REBAR_NAME);
      exit(0);
    default:
      cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
  }

}

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  // init mesh
  string filename = "/usr/local/igl/libigl/examples/shared/decimated-knight.obj" ;
  if(argc > 1)
  {
    filename = argv[1];
  }
  if(!read_triangle_mesh(filename,V,F))
  {
    return 1;
  }

  // Compute normals, centroid, colors, bounding box diagonal
  per_face_normals(V,F,N);
  normalize_row_lengths(N,N);
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  bbd =
    (V.colwise().maxCoeff() -
    V.colwise().minCoeff()).maxCoeff();

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
  rebar.load(REBAR_NAME);

  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("sorted primitives transparency");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);

  glutMainLoop();

  return 0;
}
