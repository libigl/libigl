#include <igl/OpenGL_convenience.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/normalize_row_lengths.h>
#include <igl/draw_mesh.h>
#include <igl/draw_floor.h>
#include <igl/quat_to_mat.h>
#include <igl/report_gl_error.h>
#include <igl/read.h>
#include <igl/trackball.h>
#include <igl/material_colors.h>
#include <igl/barycenter.h>
#include <igl/matlab_format.h>
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
igl::EmbreeIntersector<Eigen::MatrixXd,Eigen::MatrixXi,Eigen::Vector3d> ei;
// Running ambient occlusion
Eigen::VectorXd S;
int tot_num_samples = 0;

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

// Set up projection and model view of scene
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
    ambient_occlusion(ei,V,N,num_samples,Si);
    S *= (double)tot_num_samples;
    S += Si*(double)num_samples;
    tot_num_samples += num_samples;
    S /= (double)tot_num_samples;
    // Convert to 1-intensity
    C.resize(S.rows(),3);
    C<<S,S,S;
    C.array() = (1.0-C.array());
  }

  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // All smooth points
  glEnable( GL_POINT_SMOOTH );

  glDisable(GL_LIGHTING);
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
  draw_mesh(V,F,N,C);

  pop_object();

  // Draw a nice floor
  glPushMatrix();
  const double floor_scale = 2./bbd;
  const double floor_offset =
    -2./bbd*(V.col(1).minCoeff()+mid(1));
  glTranslated(0,floor_offset,0);
  const float GREY[4] = {0.5,0.5,0.6,1.0};
  const float DARK_GREY[4] = {0.2,0.2,0.3,1.0};
  draw_floor(GREY,DARK_GREY);
  glPopMatrix();

  pop_scene();

  report_gl_error();

  glutSwapBuffers();
  glutPostRedisplay();
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
  using namespace std;
  switch(key)
  {
    // Ctrl-c and esc exit
    case char(3):
    case char(27):
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
  if(!read("../shared/beast.obj",V,F))
  {
    return 1;
  }
  // Compute normals, centroid, colors, bounding box diagonal
  per_vertex_normals(V,F,N);
  mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
  bbd = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();

  // Init embree
  ei = EmbreeIntersector<MatrixXd,MatrixXi,Vector3d>(V,F);

  // Init glut
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("ambient-occlusion");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutMainLoop();
  return 0;
}
