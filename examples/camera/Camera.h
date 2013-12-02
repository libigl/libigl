#include <Eigen/Geometry>
#include <Eigen/Core>

class Camera
{
  public:
    //  m_zoom   Zoom of camera lens {1}
    //  m_angle  Field of view angle in degrees {45}
    //  m_aspect  Aspect ratio {1}
    //  m_near  near clipping plane {1e-2}
    //  m_far  far clipping plane {100}
    //  m_at_dist  distance of looking at point {1}
    //  m_rotation  Rotation part of rigid transformation of camera {identity}.
    //    Note that this seems to the inverse of what's stored in the old
    //    igl::Camera class.
    //  m_translation  Translation part of rigid transformation of camera
    //    {(0,0,1)}
    double m_zoom, m_angle, m_aspect, m_near, m_far, m_at_dist;
    Eigen::Quaterniond m_rotation;
    Eigen::Vector3d m_translation;
    Camera();
    // Return projection matrix that takes relative camera coordinates and
    // transforms it to viewport coordinates
    //
    // Note:
    //
    //     gluPerspective(m_angle,m_aspect,m_near,m_far);
    //
    // Is equivalent to
    //
    //     glMultMatrixd(v_camera.projection().data());
    //
    Eigen::Matrix4d projection() const;
    // Return an Affine transformation (rigid actually) that takes a world 3d coordinate and
    // transforms it into the relative camera coordinates.
    Eigen::Affine3d affine() const;
    // Return an Affine transformation (rigid actually) that takes relative
    // coordinates and tramsforms them into world 3d coordinates.
    //
    // Note:
    //
    //     gluLookAt(
    //       eye()(0), eye()(1), eye()(2),
    //       at()(0), at()(1), at()(2),
    //       up()(0), up()(1), up()(2));
    //
    // Is equivalent to
    //
    //     glMultMatrixd(camera.affine().matrix().data());
    //
    // See also: affine, eye, at, up
    Eigen::Affine3d inverse() const;
    // Returns world coordinates position of center or "eye" of camera.
    Eigen::Vector3d eye() const;
    // Returns world coordinate position of a point "eye" is looking at.
    Eigen::Vector3d at() const;
    // Returns world coordinate unit vector of "up" vector
    Eigen::Vector3d up() const;
    // Return top right corner of unit plane in relative coordinates, that is
    // (w/2,h/2,1)
    Eigen::Vector3d unit_plane() const;
    // Move dv in the relative coordinate frame of the camera (move the FPS)
    //
    // Inputs:
    //   dv  (x,y,z) displacement vector
    //
    void dolly(const Eigen::Vector3d & dv);
    // "Scale zoom": Move `eye`, but leave `at`
    //
    // Input:
    //   s  amount to scale distance to at
    void push_away(const double s);
    // Aka "Hitchcock", "Vertigo", "Spielberg" or "Trombone" zoom:
    // simultaneously dolly while changing angle so that `at` not only stays
    // put in relative coordinates but also projected coordinates. That is
    //
    // Inputs:
    //   da  change in angle in degrees
    void dolly_zoom(const double da);
    // Turn around eye so that rotation is now q
    //
    // Inputs:
    //   q  new rotation as quaternion
    void turn_eye(const Eigen::Quaterniond & q);
    // Orbit around at so that rotation is now q
    //
    // Inputs:
    //   q  new rotation as quaternion
    void orbit(const Eigen::Quaterniond & q);
    // Rotate and translate so that camera is situated at "eye" looking at "at"
    // with "up" pointing up.
    //
    // Inputs:
    //   eye  (x,y,z) coordinates of eye position
    //   at   (x,y,z) coordinates of at position
    //   up   (x,y,z) coordinates of up vector
    void look_at(
      const Eigen::Vector3d & eye,
      const Eigen::Vector3d & at,
      const Eigen::Vector3d & up);
};

// Implementation
#include <igl/PI.h>
#include <igl/EPS.h>
#include <cmath>
#include <cassert>

Camera::Camera():
  m_zoom(1),m_angle(45.0),m_aspect(1),m_near(1e-2),m_far(100),m_at_dist(1),
  m_rotation(1,0,0,0),
  m_translation(0,0,1)
{
}

Eigen::Matrix4d Camera::projection() const
{
  Eigen::Matrix4d P;
  using namespace std;
  using namespace igl;
  // http://stackoverflow.com/a/3738696/148668
  const double yScale = tan(PI*0.5 - 0.5*m_angle*PI/180.);
  // http://stackoverflow.com/a/14975139/148668
  const double xScale = yScale/m_aspect;
  P<< 
    xScale, 0, 0, 0,
    0, yScale, 0, 0,
    0, 0, -(m_far+m_near)/(m_far-m_near), -1,
    0, 0, -2.*m_near*m_far/(m_far-m_near), 0;
  return P.transpose();
}

Eigen::Affine3d Camera::affine() const
{
  using namespace Eigen;
  Affine3d t = Affine3d::Identity();
  t.rotate(m_rotation);
  t.translate(m_translation);
  return t;
}

Eigen::Affine3d Camera::inverse() const
{
  using namespace Eigen;
  Affine3d t = Affine3d::Identity();
  t.translate(-m_translation);
  t.rotate(m_rotation.conjugate());
  return t;
}

Eigen::Vector3d Camera::eye() const
{
  using namespace Eigen;
  return affine() * Vector3d(0,0,0);
}

Eigen::Vector3d Camera::at() const
{
  using namespace Eigen;
  return affine() * (Vector3d(0,0,-1)*m_at_dist);
}

Eigen::Vector3d Camera::up() const
{
  using namespace Eigen;
  Affine3d t = Affine3d::Identity();
  t.rotate(m_rotation);
  return t * Vector3d(0,1,0);
}

Eigen::Vector3d Camera::unit_plane() const
{
  using namespace igl;
  // Distance of center pixel to eye
  const double d = 1.0;
  const double a = m_aspect;
  const double theta = m_angle*PI/180.;
  const double w =
    2.*sqrt(-d*d/(a*a*pow(tan(0.5*theta),2.)-1.))*a*tan(0.5*theta);
  const double h = w/a;
  return Eigen::Vector3d(w*0.5,h*0.5,-d);
}

void Camera::dolly(const Eigen::Vector3d & dv)
{
  m_translation += dv;
}

void Camera::push_away(const double s)
{
  using namespace Eigen;
  using namespace igl;
#ifndef NDEBUG
  Vector3d old_at = at();
#endif
  const double old_at_dist = m_at_dist;
  m_at_dist = old_at_dist * s;
  dolly(Vector3d(0,0,1)*(m_at_dist - old_at_dist));
  assert((old_at-at()).squaredNorm() < DOUBLE_EPS);
}

void Camera::dolly_zoom(const double da)
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
#ifndef NDEBUG
  Vector3d old_at = at();
#endif
  const double old_angle = m_angle;
  m_angle += da;
  m_angle = min(89.,max(5.,m_angle));
  const double s = 
    (2.*tan(old_angle/2./180.*M_PI)) /
    (2.*tan(m_angle/2./180.*M_PI)) ;
  const double old_at_dist = m_at_dist;
  m_at_dist = old_at_dist * s;
  dolly(Vector3d(0,0,1)*(m_at_dist - old_at_dist));
  assert((old_at-at()).squaredNorm() < DOUBLE_EPS);
}

void Camera::turn_eye(const Eigen::Quaterniond & q)
{
  using namespace Eigen;
  using namespace igl;
  Vector3d old_eye = eye();
  // eye should be fixed
  //
  // eye_1 = R_1 * t_1 = eye_0
  // t_1 = R_1' * eye_0
  m_rotation = q;
  m_translation = m_rotation.conjugate() * old_eye;
  assert((old_eye - eye()).squaredNorm() < DOUBLE_EPS);
}

void Camera::orbit(const Eigen::Quaterniond & q)
{
  using namespace Eigen;
  using namespace igl;
#ifndef NDEBUG
  Vector3d old_at = at();
#endif
  // at should be fixed
  //
  // at_1 = R_1 * t_1 - R_1 * z = at_0
  // t_1 = R_1' * (at_0 + R_1 * z)
  m_rotation = q;
  m_translation = 
    m_rotation.conjugate() * 
      (old_at + 
         m_rotation * Vector3d(0,0,1) * m_at_dist);
  assert((old_at - at()).squaredNorm() < DOUBLE_EPS);
}

void Camera::look_at(
  const Eigen::Vector3d & eye,
  const Eigen::Vector3d & at,
  const Eigen::Vector3d & up)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  // http://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
  // Normalize vector from at to eye
  Vector3d F = eye-at;
  m_at_dist = F.norm();
  F.normalize();
  // Project up onto plane orthogonal to F and normalize
  const Vector3d proj_up = (up-(up.dot(F))*F).normalized();
  Quaterniond a,b;
  a.setFromTwoVectors(Vector3d(0,0,-1),-F);
  b.setFromTwoVectors(a*Vector3d(0,1,0),proj_up);
  m_rotation = a*b;
  m_translation = m_rotation.conjugate() * eye;
  assert(           (eye-this->eye()).squaredNorm() < DOUBLE_EPS);
  assert((F-(this->eye()-this->at()).normalized()).squaredNorm() < 
    DOUBLE_EPS);
  assert(           (at-this->at()).squaredNorm() < DOUBLE_EPS);
  assert(        (proj_up-this->up()).squaredNorm() < DOUBLE_EPS);
}
