#include <Eigen/Core>
// Draw a skeleton
//
// Inputs:
//   C  #C by dim List of joint rest positions
//   BE  #BE by 2 list of bone edge indices into C
template <typename DerivedC, typename DerivedBE>
void draw_skeleton(
  const Eigen::PlainObjectBase<DerivedC> & C,
  const Eigen::PlainObjectBase<DerivedBE> & BE);
template <typename DerivedC, typename DerivedBE, typename DerivedT>
void draw_skeleton(
  const Eigen::PlainObjectBase<DerivedC> & C,
  const Eigen::PlainObjectBase<DerivedBE> & BE,
  const Eigen::PlainObjectBase<DerivedT> & T);

// Implementation
#include <igl/PI.h>
#include <igl/OpenGL_convenience.h>
#include <Eigen/Geometry>
#include <iostream>
template <typename DerivedC, typename DerivedBE>
void draw_skeleton(
  const Eigen::PlainObjectBase<DerivedC> & C,
  const Eigen::PlainObjectBase<DerivedBE> & BE)
{
  using namespace Eigen;
  typedef Eigen::Matrix<typename DerivedC::Scalar,Dynamic,Dynamic> Mat;
  Mat I = Mat::Identity(C.cols()+1,C.cols());
  Mat T = I.replicate(BE.rows(),1);
  return draw_skeleton(C,BE,T);
}

template <typename DerivedC, typename DerivedBE, typename DerivedT>
void draw_skeleton(
  const Eigen::PlainObjectBase<DerivedC> & C,
  const Eigen::PlainObjectBase<DerivedBE> & BE,
  const Eigen::PlainObjectBase<DerivedT> & T)
{
  using namespace Eigen;
  using namespace std;
  glDisable(GL_LIGHTING);
  glColor3d(70./255.,252./255.,167./255.);
  glLineWidth(2.0);

  auto draw_sphere = [](const double r)
  {
    auto draw_circle = []()
    {
      glBegin(GL_LINE_STRIP);
      for(double th = 0;th<2.0*igl::PI;th+=(2.0*igl::PI/30.0))
      {
        glVertex3d(cos(th),sin(th),0.0);
      }
      glVertex3d(cos(0),sin(0),0.0);
      glEnd();
    };
    glPushMatrix();
    glScaled(r,r,r);
    draw_circle();
    glRotated(90.0,1.0,0.0,0.0);
    draw_circle();
    glRotated(90.0,0.0,1.0,0.0);
    draw_circle();
    glPopMatrix();
  };
  auto draw_pyramid = [](const double r)
  {
    glBegin(GL_LINE_STRIP);
    glVertex3d(0, 1,-1);
    glVertex3d(0,-1,-1);
    glVertex3d(0,-1, 1);
    glVertex3d(0, 1, 1);
    glVertex3d(0, 1,-1);
    glEnd();
    glBegin(GL_LINES);
    glVertex3d(0, 1,-1);
    glVertex3d(1,0,0);
    glVertex3d(0,-1,-1);
    glVertex3d(1,0,0);
    glVertex3d(0,-1, 1);
    glVertex3d(1,0,0);
    glVertex3d(0, 1, 1);
    glVertex3d(1,0,0);
    glEnd();
  };

  // Loop over bones
  for(int e = 0;e < BE.rows();e++)
  {
    // Draw a sphere
    auto s = C.row(BE(e,0));
    auto d = C.row(BE(e,1));
    auto b = (d-s).transpose().eval();
    double r = 0.02;
    Matrix4d Te = Matrix4d::Identity();
    Te.block(0,0,3,4) = T.block(e*4,0,4,3).transpose();
    Quaterniond q;
    q.setFromTwoVectors(Vector3d(1,0,0),b);
    glPushMatrix();
    glMultMatrixd(Te.data());
    glTranslated(s(0),s(1),s(2));
    draw_sphere(r);
    const double len = b.norm()-2.*r;
    if(len>=0)
    {
      auto u = b.normalized()*r;
      glPushMatrix();
      glTranslated(u(0),u(1),u(2));
      glMultMatrixd(Affine3d(q).matrix().data());
      glScaled(b.norm()-2.*r,r,r);
      draw_pyramid(r);
      glPopMatrix();
    }
    glTranslated(b(0),b(1),b(2));
    draw_sphere(r);
    glPopMatrix();
  }
}
