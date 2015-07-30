#include <igl/opengl/OpenGL_convenience.h>
#include <igl/opengl2/draw_floor.h>
#include <igl/opengl2/draw_mesh.h>
#include <igl/normalize_row_lengths.h>
#include <igl/per_face_normals.h>
#include <igl/quat_to_mat.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/trackball.h>
#include <igl/opengl2/unproject.h>
#include <igl/embree/EmbreeIntersector.h>

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
Eigen::MatrixXd V,N,C,mean;
// Bounding box diagonal length
double bbd;
// Faces
Eigen::MatrixXi F;
// Embree intersection structure
igl::embree::EmbreeIntersector ei;
// Hits collected
std::vector<igl::embree::Hit > hits;
// Ray information, "projection screen" corners
Eigen::Vector3f win_s,s,d,dir,NW,NE,SE,SW;
// Textures and framebuffers for "projection screen"
GLuint tex_id = 0, fbo_id = 0, dfbo_id = 0;

// Initialize textures and framebuffers. Must be called if window changes
// dimension
void init_texture()
{
  using namespace igl;
  using namespace std;
  // Set up a "render-to-texture" frame buffer and texture combo
  glDeleteTextures(1,&tex_id);
  glDeleteFramebuffersEXT(1,&fbo_id);
  glDeleteFramebuffersEXT(1,&dfbo_id);
  // http://www.opengl.org/wiki/Framebuffer_Object_Examples#Quick_example.2C_render_to_texture_.282D.29
  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //NULL means reserve texture memory, but texels are undefined
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_BGRA, GL_UNSIGNED_BYTE, NULL);
  glBindTexture(GL_TEXTURE_2D, 0);
  //-------------------------
  glGenFramebuffersEXT(1, &fbo_id);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id);
  //Attach 2D texture to this FBO
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, tex_id, 0);
  glGenRenderbuffersEXT(1, &dfbo_id);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, dfbo_id);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, width, height);
  //-------------------------
  //Attach depth buffer to FBO
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, dfbo_id);
  //-------------------------
  //Does the GPU support current FBO configuration?
  GLenum status;
  status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  switch(status)
  {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
      break;
    default:
      cout<<"error"<<endl;
  }
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void reshape(int width,int height)
{
  using namespace std;
  // Save width and height
  ::width = width;
  ::height = height;
  // Re-initialize textures and frame bufferes
  init_texture();
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

  lights();
  push_scene();
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  push_object();

  if(trackball_on)
  {
    // Draw a "laser" line
    glLineWidth(3.0);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3fv(s.data());
    glColor3f(1,0,0);
    glVertex3fv(d.data());
    glEnd();

    // Draw the start and end points used for ray
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glColor3f(1,0,0);
    glVertex3fv(s.data());
    glColor3f(0,0,1);
    glVertex3fv(d.data());
    glEnd();
  }

  // Draw the model
  glEnable(GL_LIGHTING);
  igl::opengl2::draw_mesh(V,F,N,C);

  // Draw all hits
  glBegin(GL_POINTS);
  glColor3f(0,0.2,0.2);
  for(vector<igl::embree::Hit>::iterator hit = hits.begin();
      hit != hits.end();
      hit++)
  {
    const double w0 = (1.0-hit->u-hit->v);
    const double w1 = hit->u;
    const double w2 = hit->v;
    VectorXd hitP =
      w0 * V.row(F(hit->id,0)) +
      w1 * V.row(F(hit->id,1)) +
      w2 * V.row(F(hit->id,2));
    glVertex3dv(hitP.data());
  }
  glEnd();

  pop_object();

  // Draw a nice floor
  glPushMatrix();
  glEnable(GL_LIGHTING);
  glTranslated(0,-1,0);
  igl::opengl2::draw_floor();
  glPopMatrix();

  // draw a transparent "projection screen" show model at time of hit (aka
  // mouse down)
  push_object();
  if(trackball_on)
  {
    glColor4f(0,0,0,1.0);
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glVertex3fv(SW.data());
    glVertex3fv(SE.data());
    glVertex3fv(NE.data());
    glVertex3fv(NW.data());
    glEnd();

    glDisable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tex_id);
    glColor4f(1,1,1,0.7);
    glBegin(GL_QUADS);
    glTexCoord2d(0,0);
    glVertex3fv(SW.data());
    glTexCoord2d(1,0);
    glVertex3fv(SE.data());
    glTexCoord2d(1,1);
    glVertex3fv(NE.data());
    glTexCoord2d(0,1);
    glVertex3fv(NW.data());
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
  }
  pop_object();
  pop_scene();

  // Draw a faint point over mouse
  if(!trackball_on)
  {
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(1.0,0.3,0.3,0.6);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,width,0,height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPointSize(20.0);
    glBegin(GL_POINTS);
    glVertex2fv(win_s.data());
    glEnd();
  }
  igl::opengl::report_gl_error();

  glutSwapBuffers();
  glutPostRedisplay();
}

// Initialize colors to a boring green
void init_C()
{
  C.col(0).setConstant(0.4);
  C.col(1).setConstant(0.8);
  C.col(2).setConstant(0.3);
}

void mouse_move(int mouse_x, int mouse_y)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  init_C();
  glutSetCursor(GLUT_CURSOR_CROSSHAIR);
  // Push scene and object
  push_scene();
  push_object();
  // Unproject mouse at 0 depth and some positive depth
  win_s = Vector3f(mouse_x,height-mouse_y,0);
  Vector3f win_d(mouse_x,height-mouse_y,1);
  igl::opengl2::unproject(win_s,s);
  igl::opengl2::unproject(win_d,d);
  pop_object();
  pop_scene();
  igl::opengl::report_gl_error();
  // Shoot ray at igl::opengl2::unprojected mouse in view direction
  dir = d-s;
  int num_rays_shot;
  ei.intersectRay(s,dir,hits,num_rays_shot);
  for(vector<igl::embree::Hit>::iterator hit = hits.begin();
      hit != hits.end();
      hit++)
  {
    // Change color of hit faces
    C(hit->id,0) = 1;
    C(hit->id,1) = 0.4;
    C(hit->id,2) = 0.4;
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
      glutSetCursor(GLUT_CURSOR_CROSSHAIR);
      trackball_on = false;
      hits.clear();
      init_C();
      break;
    case 0:
      // be sure this has been called recently
      mouse_move(mouse_x,mouse_y);
      // down
      glutSetCursor(GLUT_CURSOR_CYCLE);
      // collect information for trackball
      trackball_on = true;
      copy(scene_rot,scene_rot+4,down_scene_rot);
      down_mouse_x = mouse_x;
      down_mouse_y = mouse_y;
      // Collect "projection screen" locations
      push_scene();
      push_object();
      // igl::opengl2::unproject corners of window
      const double depth = 0.999;
      Vector3d win_NW(    0,height,depth);
      Vector3d win_NE(width,height,depth);
      Vector3d win_SE(width,0,depth);
      Vector3d win_SW(0,0,depth);
      igl::opengl2::unproject(win_NW,NW);
      igl::opengl2::unproject(win_NE,NE);
      igl::opengl2::unproject(win_SE,SE);
      igl::opengl2::unproject(win_SW,SW);
      // render to framebuffer
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id);
      glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, dfbo_id);
      glClearColor(0,0,0,1);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      // Render the model ---> to the framebuffer attached to the "projection
      // screen" texture
      glEnable(GL_COLOR_MATERIAL);
      glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
      glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
      igl::opengl2::draw_mesh(V,F,N,C);
      pop_object();
      pop_scene();
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
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
  string filename = "../shared/decimated-knight.obj";
  if(argc < 2)
  {
    cerr<<"Usage:"<<endl<<"    ./example input.obj"<<endl;
    cout<<endl<<"Opening default mesh..."<<endl;
  }else
  {
    // Read and prepare mesh
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

  // Init embree
  ei.init(V.cast<float>(),F.cast<int>());

  // Init glut
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH)/2.0,glutGet(GLUT_SCREEN_HEIGHT));
  glutCreateWindow("embree");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc(mouse_move);
  glutMainLoop();
  return 0;
}
