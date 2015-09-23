#include "BeachBall.h"


#include <igl/opengl2/draw_beach_ball.h>
#include <igl/quat_to_mat.h>
#include <igl/opengl2/unproject_to_zero_plane.h>
#include <igl/opengl2/unproject.h>
#include <igl/opengl2/project.h>
#include <igl/quat_mult.h>
#include <igl/quat_conjugate.h>
#include <igl/trackball.h>
#include <igl/mat_to_quat.h>

#include <OpenGL/GL.h>

#include <iostream>
#include <cmath>

BeachBall::BeachBall():
  radius(1),
  is_hover(false),
  is_down(false),
  trackball_on(false),
  down_x(-1),
  down_y(-1)
{
  r[0]=r[1]=r[2]=0;r[3]=1;
  t[0]=t[1]=t[2]=0;t[3]=0;
}

void BeachBall::pushmv() const
{
  glPushMatrix();
  glLoadIdentity();
  glMultMatrixd(mv);
}
void BeachBall::popmv() const
{
  glPopMatrix();
}

void BeachBall::push() const
{
  using namespace igl;
  glPushMatrix();
  // Translate to relative origin
  glTranslated(t[0],t[1],t[2]);
  // Rotate to relative orientation
  double mat[4*4]; quat_to_mat(r,mat); glMultMatrixd(mat);
}
void BeachBall::pop() const
{
  glPopMatrix();
}

void BeachBall::draw()
{
  using namespace std;
  using namespace igl;
  // Keep track of relative identity (modelview) matrix at draw call
  glGetDoublev(GL_MODELVIEW_MATRIX,mv);
  push();
  // Dim if hovered over
  float ld[4];
  glGetLightfv(GL_LIGHT0,GL_DIFFUSE,ld);
  if(is_hover)
  {
    float ld_hover[4];
    ld_hover[0] = 0.5*ld[0];
    ld_hover[1] = 0.5*ld[1];
    ld_hover[2] = 0.5*ld[2];
    ld_hover[3] = 0.5*ld[3];
    glLightfv(GL_LIGHT0,GL_DIFFUSE,ld_hover);
  }
  // Adjust scale only after setting origin
  glPushMatrix();
  glScaled(radius,radius,radius);
  // draw oriented glyph 
  igl::opengl2::draw_beach_ball();
  // Pop scale
  glPopMatrix();
  // Reset lighting
  glLightfv(GL_LIGHT0,GL_DIFFUSE,ld);
  // Reset relative identity
  pop();
}

bool BeachBall::in(const int x,const int y) const
{
  using namespace igl;
  using namespace std;
  pushmv();
  push();
  // Now origin is center of object
  double obj[3];
  // Check if igl::opengl2::unprojected screen point is nearby
  igl::opengl2::unproject_to_zero_plane(x,y, &obj[0], &obj[1], &obj[2]);
  bool near = (obj[0]*obj[0] + obj[1]*obj[1] + obj[2]*obj[2])<radius*radius;
  pop();
  popmv();
  return near;
}

bool BeachBall::hover(const int x,const int y)
{
  return is_hover = in(x,y);
}

bool BeachBall::down(const int x,const int y)
{
  using namespace std;
  using namespace igl;
  is_down = in(x,y);
  if(is_down)
  {
    copy(r,r+4,down_r);
    copy(t,t+4,down_t);
    down_x = x;
    down_y = y;
    trackball_on = true;
  }
  return is_down;
}

bool BeachBall::drag(const int x,const int y)
{
  using namespace igl;
  using namespace std;
  if(!is_down)
  {
    return false;
  }
  if(trackball_on)
  {
    // Pick up origin location
    pushmv();
    push();
    double origin[3];
    igl::opengl2::project(0,0,0,&origin[0],&origin[1],&origin[2]);
    pop();
    popmv();
    double rot[4];
    int VP[4];
    glGetIntegerv(GL_VIEWPORT,VP);
    trackball<double>(
      VP[2],
      VP[3],
      2,
            (down_x-origin[0]+VP[2]/2),
      VP[3]-(down_y-origin[1]+VP[3]/2),
            (     x-origin[0]+VP[2]/2),
      VP[3]-(     y-origin[1]+VP[3]/2),
      rot);
    {
      // We've computed change in rotation according to this view:
      // R = mv * r, R' = rot * (mv * r)
      // But we only want new value for r:
      // R' = mv * r'
      // mv * r' = rot * (mv * r)
      // r' = mv* * rot * sr * r
      // Convert modelview matrix to quaternion
      double sr_conj[4],scene_rot[4],t1[4],t2[4];
      mat4_to_quat(mv,scene_rot);
      //cout<<"::sr: "<<
      //  BeachBall::scene_rot[0]<<" "<<
      //  BeachBall::scene_rot[1]<<" "<<
      //  BeachBall::scene_rot[2]<<" "<<
      //  BeachBall::scene_rot[3]<<" "<<
      //  endl;
      //cout<<"  sr: "<<
      //  scene_rot[0]<<" "<<
      //  scene_rot[1]<<" "<<
      //  scene_rot[2]<<" "<<
      //  scene_rot[3]<<" "<<
      //  endl;
      // Normalize to unit quaternion (rotation only)
      double len = 
        sqrt(scene_rot[0]*scene_rot[0]+
        scene_rot[1]*scene_rot[1]+
        scene_rot[2]*scene_rot[2]+
        scene_rot[3]*scene_rot[3]);
      for(int j = 0;j<4;j++)
      {
        scene_rot[j] /= len;
      }
      // TODO: Just use Eigen::Quaterniond
      quat_conjugate(scene_rot,sr_conj);
      quat_mult<double>(sr_conj,rot,t1);
      quat_mult<double>(t1,scene_rot,t2);
      quat_mult<double>(t2,down_r,r);

    }
  }else
  {
    // We want that origin follows mouse move. First define plane we
    // igl::opengl2::projecteing screen mouse movement to as perpendicular plan passing
    // through this origin.
    pushmv();
    // down_t igl::opengl2::projected to screen to get depth value
    double p[3];
    igl::opengl2::project(down_t[0],down_t[1],down_t[2],&p[0],&p[1],&p[2]);
    // igl::opengl2::unprojected down_x,down_y with down_t depth
    double du[3];
    igl::opengl2::unproject(down_x,down_y,p[2],&du[0],&du[1],&du[2]);
    // igl::opengl2::unprojected x,y with down_t depth
    double u[3];
    igl::opengl2::unproject(x,y,p[2],&u[0], &u[1], &u[2]);
    popmv();
    // Then move this origin according to igl::opengl2::project mouse displacment
    t[0] = down_t[0] + (u[0]-du[0]);
    t[1] = down_t[1] + (u[1]-du[1]);
    t[2] = down_t[2] + (u[2]-du[2]);
  }
  return is_down;
}

bool BeachBall::right_down(const int x,const int y)
{
  using namespace std;
  using namespace igl;
  is_down = in(x,y);
  if(is_down)
  {
    copy(r,r+4,down_r);
    copy(t,t+4,down_t);
    down_x = x;
    down_y = y;
  }
  return is_down;
}

bool BeachBall::up(const int /*x*/,const int /*y*/)
{
  trackball_on = false;
  return is_down = false;
}
