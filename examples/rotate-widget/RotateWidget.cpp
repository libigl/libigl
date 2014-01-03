// Implementation
#include "RotateWidget.h"
#include <igl/OpenGL_convenience.h>
#include <igl/PI.h>
#include <igl/EPS.h>
#include <igl/ray_sphere_intersect.h>
#include <igl/project.h>
#include <igl/mat_to_quat.h>
#include <igl/trackball.h>
#include <igl/unproject.h>
#include <iostream>

Eigen::Quaterniond RotateWidget::axis_q[3] = {
  Eigen::Quaterniond(Eigen::AngleAxisd(igl::PI*0.5,Eigen::Vector3d(0,1,0))),
  Eigen::Quaterniond(Eigen::AngleAxisd(igl::PI*0.5,Eigen::Vector3d(1,0,0))),
  Eigen::Quaterniond::Identity()};

Eigen::Vector3d RotateWidget::view_direction(const int x, const int y)
{
  using namespace Eigen;
  using namespace igl;
  Vector4i vp;
  glGetIntegerv(GL_VIEWPORT,vp.data());
  const int height = vp(3);
  const Vector3d win_s(x,height-y,0), win_d(x,height-y,1);
  const Vector3d s = unproject(win_s);
  const Vector3d d = unproject(win_d);
  return d-s;
}

Eigen::Vector3d RotateWidget::view_direction(const Eigen::Vector3d & pos)
{
  using namespace Eigen;
  using namespace igl;
  const Vector3d ppos = project(pos);
  return view_direction(ppos(0),ppos(1));
}

RotateWidget::RotateWidget():
  pos(0,0,0),
  rot(Eigen::Quaterniond::Identity()),
  down_rot(rot),
  down_xy(-1,-1),drag_xy(-1,-1),
  outer_radius_on_screen(91.),
  outer_over_inner(1.13684210526),
  down_type(DOWN_TYPE_NONE), 
  selected_type(DOWN_TYPE_NONE)
{
}

bool RotateWidget::intersect(const int x, const int y, Eigen::Vector3d & hit)
{
  using namespace Eigen;
  using namespace igl;
  Vector3d view = view_direction(x,y);
  Vector4i vp;
  glGetIntegerv(GL_VIEWPORT,vp.data());
  const Vector3d ppos = project(pos);
  Vector3d uxy = unproject(Vector3d(x,vp(3)-y,ppos(2)));
  double t0,t1;
  if(!ray_sphere_intersect(uxy,view,pos,unprojected_inner_radius(),t0,t1))
  {
    return false;
  }
  hit = uxy+t0*view;
  return true;
}

double RotateWidget::unprojected_inner_radius()
{
  using namespace Eigen;
  using namespace igl;
  Vector3d off,ppos,ppos_off,pos_off;
  project(pos,ppos);
  ppos_off = ppos;
  ppos_off(0) += outer_radius_on_screen/outer_over_inner;
  unproject(ppos_off,pos_off);
  return (pos-pos_off).norm();
}

Eigen::Vector3d RotateWidget::unproject_onto(const int x, const int y)
{
  using namespace Eigen;
  using namespace igl;
  Vector4i vp;
  glGetIntegerv(GL_VIEWPORT,vp.data());
  const Vector3d ppos = project(pos);
  Vector3d uxy = unproject( Vector3d(x,vp(3)-y,ppos(2)))-pos;
  uxy /= unprojected_inner_radius()*outer_over_inner*outer_over_inner;
  return uxy;
}

bool RotateWidget::down(const int x, const int y)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  down_type = DOWN_TYPE_NONE;
  selected_type = DOWN_TYPE_NONE;
  down_xy = Vector2d(x,y);
  drag_xy = down_xy;
  down_rot = rot;
  Vector3d ppos = project(pos);
  const double r = (ppos.head(2) - down_xy).norm();
  const double thresh = 3;
  if(fabs(r - outer_radius_on_screen)<thresh)
  {
    udown = unproject_onto(x,y);
    udrag = udown;
    down_type = DOWN_TYPE_OUTLINE;
    selected_type = DOWN_TYPE_OUTLINE;
    // project mouse to same depth as pos
    return true;
  }else if(r < outer_radius_on_screen/outer_over_inner+thresh*0.5)
  {
    Vector3d hit;
    const bool is_hit = intersect(down_xy(0),down_xy(1),hit);
    if(!is_hit)
    {
      cout<<"~~~!is_hit"<<endl;
    }
    auto on_meridian = [&](
      const Vector3d & hit, 
      const Quaterniond & rot, 
      const Quaterniond m,
      Vector3d & pl_hit) -> bool
    {
      // project onto rotate plane
      pl_hit = hit-pos;
      pl_hit = (m.conjugate()*rot.conjugate()*pl_hit).eval();
      pl_hit(2) = 0;
      pl_hit = (rot*m*pl_hit).eval();
      pl_hit.normalize();
      pl_hit *= unprojected_inner_radius();
      pl_hit += pos;
      return (project(pl_hit).head(2)-project(hit).head(2)).norm()<2*thresh;
    };
    udown = (hit-pos).normalized()/outer_radius_on_screen;
    udrag = udown;
    for(int a = 0;a<3;a++)
    {
      Vector3d pl_hit;
      if(on_meridian(hit,rot,Quaterniond(axis_q[a]),pl_hit))
      {
        udown = (pl_hit-pos).normalized()/outer_radius_on_screen;
        udrag = udown;
        down_type = DownType(DOWN_TYPE_X+a);
        selected_type = down_type;
        {
          Vector3d dir3 = axis_q[a].conjugate()*down_rot.conjugate()*(hit-pos);
          dir3 = AngleAxisd(-PI*0.5,Vector3d(0,0,1))*dir3;
          dir3 = (rot*axis_q[a]*dir3).eval();
          down_dir = (project((hit+dir3).eval())-project(hit)).head(2);
          down_dir.normalize();
          // flip y because y coordinate is going to be given backwards in
          // drag()
          down_dir(1) *= -1;
        }
        return true;
      }
    }
    //assert(is_hit);
    down_type = DOWN_TYPE_TRACKBALL;
    selected_type = DOWN_TYPE_TRACKBALL;
    return true;
  }else
  {
    return false;
  }
}

bool RotateWidget::drag(const int x, const int y)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  drag_xy = Vector2d(x,y);
  switch(down_type)
  {
    case DOWN_TYPE_NONE:
      return false;
    default:
    {
      const Quaterniond & q = axis_q[down_type-DOWN_TYPE_X];
      //auto ray_plane_intersect = [&](
      //  const Vector3d & o,
      //  const Vector3d & d,
      //  const Quaterniond & rot, 
      //  const Quaterniond m) -> Vector3d
      //{
      //  const Vector3d eff_o = m.conjugate()*rot.conjugate()*o;
      //  const Vector3d eff_d = m.conjugate()*rot.conjugate()*d;
      //  const double t = eff_o(2)/eff_d(2);
      //  Vector3d p = eff_o+t*eff_d;
      //  p = (rot*m*p).eval();
      //  return p;
      //};
      //Vector3d view = view_direction(x,y);
      //Vector4i vp;
      //glGetIntegerv(GL_VIEWPORT,vp.data());
      //const Vector3d ppos = project(pos);
      //Vector3d uxy = unproject(Vector3d(x,vp(3)-y,ppos(2)));
      //udrag = ray_plane_intersect(uxy-pos,view,rot,q);
      //debug_points.clear();
      //debug_points.push_back(udrag+pos);
      // project onto rotate plane

      const double dtheta = -(drag_xy - down_xy).dot(down_dir)/
        outer_radius_on_screen/outer_over_inner*PI/2.;
      Quaterniond dq(AngleAxisd(dtheta,down_rot*q*Vector3d(0,0,1)));
      rot = dq * down_rot;
      udrag = dq * udown;
      return true;
    }
    case DOWN_TYPE_OUTLINE:
      {
        Vector3d ppos = project(pos);
        // project mouse to same depth as pos
        udrag = unproject_onto(x,y);
        const Vector2d A = down_xy - ppos.head(2);
        const Vector2d B = drag_xy - ppos.head(2);
        const double dtheta = atan2(A(0)*B(1)-A(1)*B(0),A(0)*B(0)+A(1)*B(1));
        Vector3d view = view_direction(pos).normalized();
        Quaterniond dq(AngleAxisd(dtheta,view));
        rot = dq * down_rot;
      }
      return true;
    case DOWN_TYPE_TRACKBALL:
      {
        Vector3d ppos = project(pos);
        const int w = outer_radius_on_screen/outer_over_inner;//vp(2);
        const int h = w;//vp(3);
        Quaterniond dq;
        trackball(
          w,h,
          1,
          Quaterniond::Identity(),
          (down_xy(0)-ppos(0)+w/2),
          (down_xy(1)-ppos(1)+h/2),
          (     x-ppos(0)+w/2),
          (     y-ppos(1)+h/2),
          dq);
        // We've computed change in rotation according to this view:
        // R = mv * r, R' = rot * (mv * r)
        // But we only want new value for r:
        // R' = mv * r'
        // mv * r' = rot * (mv * r)
        // r' = mv* * rot * mv * r
        Matrix4d mv;
        glGetDoublev(GL_MODELVIEW_MATRIX,mv.data());
        Quaterniond scene_rot;
        // Convert modelview matrix to quaternion
        mat4_to_quat(mv.data(),scene_rot.coeffs().data());
        scene_rot.normalize();
        rot = scene_rot.conjugate() * dq * scene_rot * down_rot;
      }
      return true;
  }
}

bool RotateWidget::up(const int x, const int y)
{
  down_type = DOWN_TYPE_NONE;
  return false;
}

bool RotateWidget::is_down()
{
  return down_type != DOWN_TYPE_NONE;
}

void RotateWidget::draw()
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  const Vector4d 
    MAYA_GREEN(128./255.,242./255.,0./255.,1.),
    MAYA_YELLOW(255./255.,247./255.,50./255.,1.),
    MAYA_RED(234./255.,63./255.,52./255.,1.),
    MAYA_BLUE(0./255.,73./255.,252./255.,1.),
    MAYA_PURPLE(180./255.,73./255.,200./255.,1.),
    MAYA_GREY(0.5,0.5,0.5,1.0),
    MAYA_CYAN(131./255.,219./255.,252./255.,1.);

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glLineWidth(2.0);

  double r = unprojected_inner_radius();
  Vector3d view = view_direction(pos).normalized();

  auto draw_circle = [&](const bool cull)
  {
    Vector3d view = view_direction(pos).normalized();
    glBegin(GL_LINES);
    const double th_step = (2.0*igl::PI/100.0);
    for(double th = 0;th<2.0*igl::PI+th_step;th+=th_step)
    {
      Vector3d a(cos(th),sin(th),0.0);
      Vector3d b(cos(th+th_step),sin(th+th_step),0.0);
      if(!cull || (0.5*(a+b)).dot(view)<FLOAT_EPS)
      {
        glVertex3dv(a.data());
        glVertex3dv(b.data());
      }
    }
    glEnd();
  };


  glPointSize(5.);
  glColor4dv(MAYA_PURPLE.data());
  glBegin(GL_POINTS);
  for(auto & p : debug_points)
  {
    glVertex3dv(p.data());
  }
  glEnd();

  glPushMatrix();
  glTranslated(pos(0),pos(1),pos(2));

  glScaled(r,r,r);
  // Draw outlines
  {
    glPushMatrix();
    glColor4dv(MAYA_GREY.data());
    Quaterniond q;
    q.setFromTwoVectors(Vector3d(0,0,1),view);
    glMultMatrixd(Affine3d(q).matrix().data());
    draw_circle(false);
    glScaled(outer_over_inner,outer_over_inner,outer_over_inner);
    if(selected_type == DOWN_TYPE_OUTLINE)
    {
      glColor4dv(MAYA_YELLOW.data());
    }else
    {
      glColor4dv(MAYA_CYAN.data());
    }
    draw_circle(false);
    glPopMatrix();
  }
  // Draw quartiles
  {
    glPushMatrix();
    glMultMatrixd(Affine3d(rot).matrix().data());
    if(selected_type == DOWN_TYPE_Z)
    {
      glColor4dv(MAYA_YELLOW.data());
    }else
    {
      glColor4dv(MAYA_BLUE.data());
    }
    draw_circle(true);
    if(selected_type == DOWN_TYPE_Y)
    {
      glColor4dv(MAYA_YELLOW.data());
    }else
    {
      glColor4dv(MAYA_GREEN.data());
    }
    glRotated(90.0,1.0,0.0,0.0);
    draw_circle(true);
    if(selected_type == DOWN_TYPE_X)
    {
      glColor4dv(MAYA_YELLOW.data());
    }else
    {
      glColor4dv(MAYA_RED.data());
    }
    glRotated(90.0,0.0,1.0,0.0);
    draw_circle(true);
    glPopMatrix();
  }
  glColor3dv(MAYA_GREY.data());
  draw_guide();
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
};

void RotateWidget::draw_guide()
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  // http://www.codeproject.com/Articles/23444/A-Simple-OpenGL-Stipple-Polygon-Example-EP_OpenGL_
  const GLubyte halftone[] = {
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
    0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55};


  switch(down_type)
  {
    case DOWN_TYPE_NONE:
    case DOWN_TYPE_TRACKBALL:
      return;
    case DOWN_TYPE_OUTLINE:
      glScaled(outer_over_inner,outer_over_inner,outer_over_inner);
      break;
    default:
      break;
  }
  const Vector3d nudown(udown.normalized()), 
    nudrag(udrag.normalized());
  glPushMatrix();
  glDisable(GL_POINT_SMOOTH);
  glPointSize(5.);
  glBegin(GL_POINTS);
  glVertex3dv(nudown.data());
  glVertex3d(0,0,0);
  glVertex3dv(nudrag.data());
  glEnd();
  glBegin(GL_LINE_STRIP);
  glVertex3dv(nudown.data());
  glVertex3d(0,0,0);
  glVertex3dv(nudrag.data());
  glEnd();
  glEnable(GL_POLYGON_STIPPLE);
  glPolygonStipple(halftone);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3d(0,0,0);
  Quaterniond dq = rot * down_rot.conjugate();
  //dq.setFromTwoVectors(nudown,nudrag);
  for(double t = 0;t<1;t+=0.1)
  {
    const Vector3d p = Quaterniond::Identity().slerp(t,dq) * nudown;
    glVertex3dv(p.data());
  }
  glVertex3dv(nudrag.data());
  glEnd();
  glDisable(GL_POLYGON_STIPPLE);
  glPopMatrix();
}

