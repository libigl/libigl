#include "BeachBall.h"

#include <igl/quat_to_mat.h>
#include <igl/opengl/report_gl_error.h>
#include <igl/trackball.h>
#include <igl/canonical_quaternions.h>
#include <igl/PI.h>

#include <GLUT/GLUT.h>

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

double back[3] = {85.0/255.0,85.0/255.0,85.0/255.0};
double width,height;
double scene_rot[4] = {0,0,0,1};
double * BeachBall::scene_rot = NULL;
double down_scene_rot[4] = {0,0,0,1};
int down_mouse_x=-1,down_mouse_y=-1;
bool trackball_on = false;
const float light_pos[4] = {-0.2, -0.1, 1,0};

std::vector<BeachBall> balls;
float plane[4] = {0,0,1,1.5};

int display_count = 0;

void push_scene()
{
  using namespace igl;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //gluOrtho2D(0,width,0,height);
  gluPerspective(15,(double)width/(double)height,1e-2,300);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glPushMatrix();
  gluLookAt(0,0,20,0,0,0,0,1,0);
  double mat[4*4];
  quat_to_mat(scene_rot,mat);
  glMultMatrixd(mat);
}

void lights()
{
  using namespace std;
  using namespace igl;

  glEnable(GL_LIGHTING);
  float ones[4] = {1.0,1.0,1.0,1.0};
  glLightfv(GL_LIGHT0,GL_AMBIENT,ones);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,ones);
  glLightfv(GL_LIGHT0,GL_SPECULAR,ones);
  glLightfv(GL_LIGHT0,GL_POSITION,light_pos);
  glEnable(GL_LIGHT0);
}

void pop_scene()
{
  glPopMatrix();
}

void display()
{
  using namespace std;
  using namespace igl;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  if(true || display_count==0 || display_count==200)
  {
    glClearColor(back[0],back[1],back[2],0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    
    lights();
    push_scene();

    //// Normal list
    //for(int i = 0;i<(int)balls.size();i++)
    //{
    //  balls[i].draw();
    //}

    // single root
    balls[0].draw();
    balls[0].push();
    for(int i = 1;i<(int)balls.size();i++)
    {
      balls[i].radius = 0.5;
      balls[i].draw();
    }
    balls[0].pop();

  // Draw floor
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(10.0,1);
    float size = 10;
    glBegin(GL_QUADS);
    glNormal3f(0,0,1);
    glColor4f(0,0.2,0.8,0.5);
    glVertex4f(-size,-size,-plane[3],1);
    glColor4f(0,0,1,0.5);
    glVertex4f( size,-size,-plane[3],1);
    glColor4f(0.0,0.1,0.5,0.5);
    glVertex4f( size, size,-plane[3],1);
    //glColor4f(back[0],back[1],back[2],1);
    glColor4f(0.0,0.1,0.5,0.5);
    glVertex4f(-size, size,-plane[3],1);
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
    pop_scene();
  }


  //// Animation
  //for(int i = 0;i<(int)balls.size();i++)
  //{
  //  balls[i].r[0] += cos(double(i*display_count)/1000.0) + sin(double(display_count)/1000.0);
  //  balls[i].r[1] += sin(double(2*i*display_count)/1000.0) + cos(double(2*display_count)/1000.0);
  //  balls[i].r[2] += cos(double(3*i*display_count)/1000.0) + sin(double(3*display_count)/1000.0);
  //  balls[i].r[3] += sin(double(4*i*display_count)/1000.0) + cos(double(4*display_count)/1000.0);
  //  double len = 
  //    sqrt(balls[i].r[0]*balls[i].r[0]+
  //    balls[i].r[1]*balls[i].r[1]+
  //    balls[i].r[2]*balls[i].r[2]+
  //    balls[i].r[3]*balls[i].r[3]);
  //  for(int j = 0;j<4;j++)
  //  {
  //    balls[i].r[j] /= len;
  //  }
  //}

  igl::opengl::report_gl_error();
  glutSwapBuffers();
  glutPostRedisplay();

  display_count++;
}

void reshape(int width, int height)
{
  using namespace std;
  ::width = width;
  ::height = height;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,width,height);
}

void key(unsigned char key, int /*mouse_x*/, int /*mouse_y*/)
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

void mouse(int glutButton, int glutState, int mouse_x, int mouse_y)
{
  mouse_y = height-mouse_y;
  using namespace std;
  push_scene();
  if(glutState==1)
  {
    for(int i = 0;i<(int)balls.size();i++)
    {
      balls[i].up(mouse_x,mouse_y);
      balls[i].hover(mouse_x,mouse_y);
    }
    trackball_on = false;
  }else if(glutState==0)
  {
    // down
    bool any = false;
    for(int i = 0;!any && i<(int)balls.size();i++)
    {
      if(glutButton==0)
      {
        any |= balls[i].down(mouse_x,mouse_y);
      }else if(glutButton==2)
      {
        any |= balls[i].right_down(mouse_x,mouse_y);
      }
    }
    trackball_on = !any;
    copy(scene_rot,scene_rot+4,down_scene_rot);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
  }
  pop_scene();
}

void mouse_drag(int mouse_x, int mouse_y)
{
  mouse_y = height-mouse_y;
  using namespace std;
  using namespace igl;
  if(trackball_on)
  {
    trackball<double>(
      width,
      height,
      2,
      down_scene_rot,
      down_mouse_x,
      height-down_mouse_y,
      mouse_x,
      height-mouse_y,
      scene_rot);
  }
  push_scene();
  for(int i = 0;i<(int)balls.size();i++)
  {
    balls[i].drag(mouse_x,mouse_y);
  }
  pop_scene();

}

void mouse_move(int mouse_x, int mouse_y)
{
  mouse_y = height-mouse_y;
  push_scene();
  bool any = false;
  for(int i = 0;!any && i<(int)balls.size();i++)
  {
    any |= balls[i].hover(mouse_x,mouse_y);
  }
  pop_scene();
}


int main(int argc,char * argv[])
{
  using namespace igl;
  using namespace std;
  glutInit(&argc,argv);
  glutInitDisplayString( "rgba depth double samples>=8 ");
  glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH),glutGet(GLUT_SCREEN_HEIGHT)/2);
  glutCreateWindow(__FILE__);
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mouse);
  glutMotionFunc(mouse_drag);
  glutPassiveMotionFunc(mouse_move);

  //copy(XZ_PLANE_QUAT_D,XZ_PLANE_QUAT_D+4,scene_rot);

  balls.push_back(BeachBall());
  const int children = 10;
  for(int i = 1;i<=children;i++)
  {
    balls.push_back(BeachBall());
    balls[i].t[0] = 2.5*sin(double(i)*2.0*PI/double(children));
    balls[i].t[1] = 2.5*cos(double(i)*2.0*PI/double(children));
  }
  BeachBall::scene_rot = scene_rot;

  glutMainLoop();
  return 0;
}
