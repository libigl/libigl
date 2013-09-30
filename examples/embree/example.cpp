#include <igl/OpenGL_convenience.h>
#include <igl/per_face_normals.h>
#include <igl/read.h>
#include <igl/normalize_row_lengths.h>
#include <igl/draw_mesh.h>
#include <igl/draw_floor.h>
#include <igl/unproject.h>
#include <igl/quat_to_mat.h>
#include <igl/embree/EmbreeIntersector.h>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <Eigen/Core>

#include <iostream>

int width,height;
float scene_rot[4] = {0,0,0,1};
float light_pos[4] = {0.1,0.1,-0.9,0};
Eigen::MatrixXd V,N,C,mean;
double bbd;
Eigen::MatrixXi F;
igl::EmbreeIntersector<Eigen::MatrixXd,Eigen::MatrixXi,Eigen::Vector3d> ei;
// Ray
Eigen::Vector3d s,d,dir;

void reshape(int width,int height)
{
  using namespace std;
  ::width = width;
  ::height = height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,width,height);
}

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

void push_scene()
{
  using namespace igl;
  //gluOrtho2D(0,width,0,height);
  gluPerspective(45,(double)width/(double)height,1e-2,100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0,0,3,0,0,0,0,1,0);
  glPushMatrix();
  float mat[4*4];
  quat_to_mat(scene_rot,mat);
  glMultMatrixf(mat);
}

void push_object()
{
  glPushMatrix();
  glScaled(2./bbd,2./bbd,2./bbd);
  glTranslated(-mean(0,0),-mean(0,1),-mean(0,2));
}

void pop_scene()
{
  glPopMatrix();
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
  glClearColor(back[0],back[1],back[2],0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  lights();
  push_scene();

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  //glColorMaterial(GL_FRONT, GL_DIFFUSE);
  //glColorMaterial(GL_FRONT, GL_AMBIENT);
  //glColorMaterial(GL_FRONT, GL_SPECULAR);


  push_object();
  draw_mesh(V,F,N,C);

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  glColor3f(1,0,0);
  glVertex3dv(s.data());
  glColor3f(0,0,1);
  glVertex3dv(d.data());
  glEnd();
  Vector3d n,f;
  n = s+1000.0*dir;
  f = d-1000.0*dir;
  glBegin(GL_LINE);
  glColor3f(1,0,0);
  glVertex3dv(n.data());
  glColor3f(1,0,0);
  glVertex3dv(f.data());
  glEnd();

  pop_object();

  glPushMatrix();
  glEnable(GL_LIGHTING);
  glTranslated(0,-1,0);
  draw_floor();
  glPopMatrix();

  pop_scene();

  glutSwapBuffers();
  glutPostRedisplay();
}

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
  init_C();
  push_scene();
  push_object();
  Vector3d win_s(mouse_x,height-mouse_y,0);
  Vector3d win_d(mouse_x,height-mouse_y,1);
  unproject(win_s,s);
  unproject(win_d,d);
  dir = d-s;
  embree::Hit hit;
  if(ei.intersectRay(s,d,hit))
  {
    cout<<"hit!"<<endl;
    C(hit.id0 % F.rows(),0) = 1;
    C(hit.id0 % F.rows(),1) = 0.4;
    C(hit.id0 % F.rows(),2) = 0.4;
  }
  pop_object();
  pop_scene();
}


void key(unsigned char key, int mouse_x, int mouse_y)
{
  using namespace std;
  switch(key)
  {
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
  if(!read("../shared/cheburashka.obj",V,F))
  {
    return 1;
  }
  per_face_normals(V,F,N);
  mean = V.colwise().mean();
  C.resize(F.rows(),3);
  init_C();
  bbd = 
    (V.colwise().maxCoeff() -
    V.colwise().minCoeff()).maxCoeff();
  normalize_row_lengths(N,N);

  // Init embree
  cout<<"Flipping faces..."<<endl;
  MatrixXi FF;
  FF.resize(F.rows()*2,F.cols());
  FF << F, F.rowwise().reverse().eval();
  cout<<"Initializing Embree..."<<endl;
  ei = EmbreeIntersector<MatrixXd,MatrixXi,Vector3d>(V,FF);

  // Init glut
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH),glutGet(GLUT_SCREEN_HEIGHT));
  glutInitWindowSize(450,300);
  glutCreateWindow("embree");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutPassiveMotionFunc(mouse_move);
  glutMainLoop();
  return 0;
}
