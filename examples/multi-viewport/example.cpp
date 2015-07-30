#include <igl/Camera.h>
#include <igl/EPS.h>
#include <igl/opengl/OpenGL_convenience.h>
#include <igl/STR.h>
#include <igl/Viewport.h>
#include <igl/canonical_quaternions.h>
#include <igl/opengl2/draw_beach_ball.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/normalize_row_lengths.h>
#include <igl/per_face_normals.h>
#include <igl/opengl2/project.h>
#include <igl/quat_to_mat.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/trackball.h>
#include <igl/opengl2/unproject.h>
#include <igl/opengl2/unproject_to_zero_plane.h>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <Eigen/Core>

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#define IGL_HEADER_ONLY
#include <igl/opengl2/draw_floor.h>

#define NUM_VIEWPORTS 4
class AugViewport : public igl::Viewport
{
  public:
    igl::Camera camera;
} viewports[NUM_VIEWPORTS];
double horiz = 0.5;
double vert = 0.5;
bool horiz_on = false, vert_on = false;

// Width and height of window
int width,height;
// information at mouse down
igl::Camera down_camera;
bool trackball_on = false;
int down_mouse_x,down_mouse_y,move_x,move_y;
int down_vp;
// Position of light
float light_pos[4] = {0.1,0.1,-0.9,0};
// Vertex positions, normals, colors and centroid
Eigen::MatrixXd V,N,C,mean;
// Bounding box diagonal length
double bbd;
// Faces
Eigen::MatrixXi F;
Eigen::Vector3d ball;


void init_viewports()
{
  using namespace igl;
  using namespace std;
  for(auto & vp : viewports)
  {
    vp.camera.push_away(5.);
  }
  viewports[0].camera.dolly_zoom(0.-viewports[0].camera.m_angle);
  viewports[1].camera.dolly_zoom(0.-viewports[1].camera.m_angle);
  viewports[2].camera.dolly_zoom(0.-viewports[2].camera.m_angle);
  viewports[3].camera.dolly_zoom(25.-viewports[3].camera.m_angle);
  // Above view
  double XZ_PLANE_QUAT_D_FLIP[4];
  copy(XZ_PLANE_QUAT_D,XZ_PLANE_QUAT_D+4,XZ_PLANE_QUAT_D_FLIP);
  XZ_PLANE_QUAT_D_FLIP[0] *= -1.0;
  // Straight on
  copy(
    XZ_PLANE_QUAT_D_FLIP,
    XZ_PLANE_QUAT_D_FLIP+4,
    viewports[0].camera.m_rotation_conj.coeffs().data());
  // Left side view
  copy(
    XY_PLANE_QUAT_D,
    XY_PLANE_QUAT_D+4,
    viewports[1].camera.m_rotation_conj.coeffs().data());
  copy(
    CANONICAL_VIEW_QUAT_D[14],
    CANONICAL_VIEW_QUAT_D[14]+4,
    viewports[2].camera.m_rotation_conj.coeffs().data());
  // Straight on
  copy(
    XY_PLANE_QUAT_D,
    XY_PLANE_QUAT_D+4,
    viewports[3].camera.m_rotation_conj.coeffs().data());
}

const double BAR_THICKNESS = 3.0;
void clamp(const double horiz, const double vert, double & eff_h, double & eff_v)
{
  eff_h = horiz;
  eff_v = vert;
  const double MIN_H = (BAR_THICKNESS/2.0)/(double)height;
  const double MIN_V = (BAR_THICKNESS/2.0)/(double)width;
  eff_h = eff_h < MIN_H ? MIN_H : eff_h;
  eff_v = eff_v < MIN_V ? MIN_V : eff_v;
  eff_h = eff_h > (1.0-MIN_H) ? (1.0-MIN_H) : eff_h;
  eff_v = eff_v > (1.0-MIN_V) ? (1.0-MIN_V) : eff_v;
}

// Viewports are arranged like planar quadrants (CCW)
// /-----.-----\
// |  1  |  0  |
// -------------
// |  2  |  3  |
// \-----.-----/
void reshape_viewports()
{
  // Make horiz and vert sane
  horiz = (horiz < 0 ? 0 : horiz);
  vert = (vert < 0 ? 0 : vert);
  horiz = (horiz > 1 ? 1 : horiz);
  vert = (vert > 1 ? 1 : vert);
  // Make effective horiz and vert even saner
  double eff_h,eff_v;
  clamp(horiz,vert,eff_h,eff_v);
  viewports[0].reshape(eff_v*width,eff_h*height,(1.-eff_v)*width,(1.-eff_h)*height);
  viewports[1].reshape(         0,eff_h*height,     eff_v*width,(1.-eff_h)*height);
  viewports[2].reshape(         0,           0,     eff_v*width,eff_h*height);
  viewports[3].reshape(eff_v*width,           0,(1.-eff_v)*width,eff_h*height);
  for(auto & vp : viewports)
  {
    vp.camera.m_aspect = (double)vp.width/(double)vp.height;
  }
}

void reshape(int width,int height)
{
  using namespace std;
  // Save width and height
  ::width = width;
  ::height = height;
  reshape_viewports();
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
void push_scene(const AugViewport & vp)
{
  using namespace igl;
  using namespace std;
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  auto & camera = vp.camera;
  glMultMatrixd(camera.projection().data());
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

// Scale and shift for object
void push_object()
{
  glPushMatrix();
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-mean(0,0),-mean(0,1),-mean(0,2));
}

void pop_object()
{
  glPopMatrix();
}

const float back[4] = {190.0/255.0,190.0/255.0,190.0/255.0,0};
void display()
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // All smooth points
  glEnable( GL_POINT_SMOOTH );

  // Flashlights
  lights();
  for(int vp = 0;vp<NUM_VIEWPORTS;vp++)
  {
    if(
      viewports[vp].width <= 0 ||
      viewports[vp].height <= 0)
    {
      continue;
    }
    glViewport(
      viewports[vp].x,
      viewports[vp].y,
      viewports[vp].width,
      viewports[vp].height);
    push_scene(viewports[vp]);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    push_object();
    // Draw the model
    glEnable(GL_LIGHTING);
    igl::opengl2::draw_mesh(V,F,N,C);
    pop_object();

    // Draw a nice floor
    glPushMatrix();
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glTranslated(0,-1,0);
    if(igl::opengl2::project(Vector3d(0,0,0))(2) - igl::opengl2::project(Vector3d(0,1,0))(2) > -FLOAT_EPS)
    {
      igl::opengl2::draw_floor_outline();
    }
    igl::opengl2::draw_floor();
    glPopMatrix();
    glDisable(GL_CULL_FACE);

    // Draw a ball at the mouse
    if(!horiz_on && !vert_on && !trackball_on)
    {
      glPushMatrix();
      glTranslated(ball(0),ball(1),ball(2));
      glScaled(0.1,0.1,0.1);
      igl::opengl2::draw_beach_ball();
      glPopMatrix();
    }

    pop_scene();
  }

  // Screen space
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glViewport(0,0,width,height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,width,0,height);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Display mouse position at cursor
  string str;
  void * font = GLUT_BITMAP_HELVETICA_18;
  glColor3f(1,0,0);
  glRasterPos2d(move_x+18,height-move_y+18);
  str = STR("("<<move_x<<","<<move_y<<")");
  for_each(str.begin(),str.end(),bind1st(ptr_fun(&glutBitmapCharacter),font));
  glColor3f(0,0,1);
  glRasterPos2d(move_x+18,height-move_y-1.5*18);
  for_each(str.begin(),str.end(),bind1st(ptr_fun(&glutBitmapCharacter),font));

  double eff_h,eff_v;
  clamp(horiz,vert,eff_h,eff_v);
  glLineWidth(BAR_THICKNESS);
  glColor3f(0.5,0.5,0.5);
  glBegin(GL_LINES);
  glVertex2f(0,eff_h*height);
  glVertex2f(width,eff_h*height);
  glVertex2f(eff_v*width,0);
  glVertex2f(eff_v*width,height);
  glEnd();

  glLineWidth(BAR_THICKNESS/3.0);
  glColor3f(0.8,0.8,0.8);
  glBegin(GL_LINES);
  glVertex2f(0,eff_h*height);
  glVertex2f(width,eff_h*height);
  glVertex2f(eff_v*width,0);
  glVertex2f(eff_v*width,height);
  glEnd();


  igl::opengl::report_gl_error();

  glutSwapBuffers();
}

// Initialize colors to a boring green
void init_C()
{
  C.col(0).setConstant(0.4);
  C.col(1).setConstant(0.8);
  C.col(2).setConstant(0.3);
}

const double SNAP_DIST = 10;
void mouse_move(int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  move_x = mouse_x;
  move_y = mouse_y;
  const int in_vp = find_if(viewports,viewports+NUM_VIEWPORTS,
    [&mouse_x,&mouse_y](const AugViewport & vp) -> bool
    {return vp.inside(mouse_x,height-mouse_y);})-viewports;
  if(in_vp < NUM_VIEWPORTS)
  {
    if(
      viewports[in_vp].width > 0 &&
      viewports[in_vp].height > 0)
    {
      glViewport(
        viewports[in_vp].x,
        viewports[in_vp].y,
        viewports[in_vp].width,
        viewports[in_vp].height);
      push_scene(viewports[in_vp]);
      Vector3d screen_ball(mouse_x,height-mouse_y,0);
      igl::opengl2::unproject_to_zero_plane(screen_ball,ball);
      pop_scene();
    }
  }
  if( (fabs((height-mouse_y) - horiz*height) < 2.*SNAP_DIST)
   && (fabs(mouse_x - vert*width) < 2.*SNAP_DIST))
  {
    glutSetCursor(GLUT_CURSOR_TOP_LEFT_CORNER);
  } else if( fabs((height-mouse_y) - horiz*height) < SNAP_DIST)
  {
    glutSetCursor(GLUT_CURSOR_UP_DOWN);
  } else if( fabs(mouse_x - vert*width) < SNAP_DIST)
  {
    glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
  } else
  {
    glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
  }
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
      horiz_on = vert_on = false;
      if( (fabs((height-mouse_y) - horiz*height) < 2.*SNAP_DIST)
          && (fabs(mouse_x - vert*width) < 2.*SNAP_DIST))
      {
        glutSetCursor(GLUT_CURSOR_TOP_LEFT_CORNER);
        horiz_on = vert_on = true;
      } else if( fabs((height-mouse_y) - horiz*height) < SNAP_DIST)
      {
        glutSetCursor(GLUT_CURSOR_UP_DOWN);
        horiz_on = true;
      } else if( fabs(mouse_x - vert*width) < SNAP_DIST)
      {
        glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
        vert_on = true;
      } else
      {
        down_vp = find_if(viewports,viewports+NUM_VIEWPORTS,
          [&mouse_x,&mouse_y](const AugViewport & vp) -> bool
          {return vp.inside(mouse_x,height-mouse_y);})-viewports;
        // down
        if(down_vp < NUM_VIEWPORTS)
        {
          glutSetCursor(GLUT_CURSOR_CYCLE);
          // collect information for trackball
          trackball_on = true;
          down_camera = viewports[down_vp].camera;
          down_mouse_x = mouse_x;
          down_mouse_y = mouse_y;
        }
      }
    break;
  }
  glutPostRedisplay();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  using namespace igl;
  mouse_move(mouse_x,mouse_y);

  if(horiz_on)
  {
    glutSetCursor(GLUT_CURSOR_UP_DOWN);
    horiz = (double)(height-mouse_y)/(double)height;
    reshape_viewports();
  }
  if(vert_on)
  {
    if(horiz_on)
    {
      glutSetCursor(GLUT_CURSOR_TOP_LEFT_CORNER);
    }else
    {
      glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
    }
    vert = (double)mouse_x/(double)width;
    reshape_viewports();
  }
  if(trackball_on)
  {
    glutSetCursor(GLUT_CURSOR_CYCLE);
    // Rotate according to trackball
    trackball<double>(
      viewports[down_vp].width,
      viewports[down_vp].height,
      2.0,
      down_camera.m_rotation_conj.coeffs().data(),
      viewports[down_vp].mouse_x(down_mouse_x),
      viewports[down_vp].mouse_y(down_mouse_y,height),
      viewports[down_vp].mouse_x(mouse_x),
      viewports[down_vp].mouse_y(mouse_y,height),
      viewports[down_vp].camera.m_rotation_conj.coeffs().data());
  }
  glutPostRedisplay();
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
    case 'Z':
    {
      const int in_vp = find_if(viewports,viewports+NUM_VIEWPORTS,
          [&mouse_x,&mouse_y](const AugViewport & vp) -> bool
          {return vp.inside(mouse_x,height-mouse_y);})-viewports;
      if(in_vp < NUM_VIEWPORTS)
      {
        igl::snap_to_canonical_view_quat(
          viewports[in_vp].camera.m_rotation_conj.coeffs().data(),
          1.0,
          viewports[in_vp].camera.m_rotation_conj.coeffs().data());
      }
      break;
    }
    default:
      cout<<"Unknown key command: "<<key<<" "<<int(key)<<endl;
  }

  glutPostRedisplay();
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
  mean = V.colwise().mean();
  C.resize(F.rows(),3);
  init_C();
  bbd =
    (V.colwise().maxCoeff() -
    V.colwise().minCoeff()).maxCoeff();


  // Init viewports
  init_viewports();

  // Init glut
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("multi-viewport");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc(mouse_move);
  glutMainLoop();
  return 0;
}
