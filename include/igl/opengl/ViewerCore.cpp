// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerCore.h"
#include "gl.h"
#include "../quat_to_mat.h"
#include "../snap_to_fixed_up.h"
#include "../look_at.h"
#include "../frustum.h"
#include "../ortho.h"
#include "../massmatrix.h"
#include "../barycenter.h"
#include "../PI.h"
#include <Eigen/Geometry>
#include <iostream>

IGL_INLINE void igl::opengl::ViewerCore::align_camera_center(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F)
{
  if(V.rows() == 0)
    return;

  get_scale_and_shift_to_fit_mesh(V,F,model_zoom,model_translation);
  // Rather than crash on empty mesh...
  if(V.size() > 0)
  {
    object_scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  }
}

IGL_INLINE void igl::opengl::ViewerCore::get_scale_and_shift_to_fit_mesh(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  float& zoom,
  Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  Eigen::MatrixXd BC;
  if (F.rows() <= 1)
  {
    BC = V;
  } else
  {
    igl::barycenter(V,F,BC);
  }
  return get_scale_and_shift_to_fit_mesh(BC,zoom,shift);
}

IGL_INLINE void igl::opengl::ViewerCore::align_camera_center(
  const Eigen::MatrixXd& V)
{
  if(V.rows() == 0)
    return;

  get_scale_and_shift_to_fit_mesh(V,model_zoom,model_translation);
  // Rather than crash on empty mesh...
  if(V.size() > 0)
  {
    object_scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  }
}

IGL_INLINE void igl::opengl::ViewerCore::get_scale_and_shift_to_fit_mesh(
  const Eigen::MatrixXd& V,
  float& zoom,
  Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  auto min_point = V.colwise().minCoeff();
  auto max_point = V.colwise().maxCoeff();
  auto centroid  = (0.5*(min_point + max_point)).eval();
  shift.setConstant(0);
  shift.head(centroid.size()) = -centroid.cast<float>();
  zoom = 2.0 / (max_point-min_point).array().abs().maxCoeff();
}


IGL_INLINE void igl::opengl::ViewerCore::clear_framebuffers()
{
  glClearColor(background_color[0],
               background_color[1],
               background_color[2],
               1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

IGL_INLINE void igl::opengl::ViewerCore::draw(
  ViewerData& data,
  bool update_matrices)
{
  using namespace std;
  using namespace Eigen;

  if (depth_test)
    glEnable(GL_DEPTH_TEST);
  else
    glDisable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  /* Bind and potentially refresh mesh/line/point data */
  if (data.dirty)
  {
    data.updateGL(data, data.invert_normals,data.meshgl);
    data.dirty = MeshGL::DIRTY_NONE;
  }
  data.meshgl.bind_mesh();

  // Initialize uniform
  glViewport(viewport(0), viewport(1), viewport(2), viewport(3));

  if(update_matrices)
  {
    model = Eigen::Matrix4f::Identity();
    view  = Eigen::Matrix4f::Identity();
    proj  = Eigen::Matrix4f::Identity();

    // Set view
    look_at( camera_eye, camera_center, camera_up, view);

    float width  = viewport(2);
    float height = viewport(3);

    // Set projection
    if (orthographic)
    {
      float length = (camera_eye - camera_center).norm();
      float h = tan(camera_view_angle/360.0 * igl::PI) * (length);
      ortho(-h*width/height, h*width/height, -h, h, camera_dnear, camera_dfar,proj);
    }
    else
    {
      float fH = tan(camera_view_angle / 360.0 * igl::PI) * camera_dnear;
      float fW = fH * (double)width/(double)height;
      frustum(-fW, fW, -fH, fH, camera_dnear, camera_dfar,proj);
    }
    // end projection

    // Set model transformation
    float mat[16];
    igl::quat_to_mat(trackball_angle.coeffs().data(), mat);

    for (unsigned i=0;i<4;++i)
      for (unsigned j=0;j<4;++j)
        model(i,j) = mat[i+4*j];

    // Why not just use Eigen::Transform<double,3,Projective> for model...?
    model.topLeftCorner(3,3)*=camera_zoom;
    model.topLeftCorner(3,3)*=model_zoom;
    model.col(3).head(3) += model.topLeftCorner(3,3)*model_translation;
  }

  // Send transformations to the GPU
  GLint modeli = glGetUniformLocation(data.meshgl.shader_mesh,"model");
  GLint viewi  = glGetUniformLocation(data.meshgl.shader_mesh,"view");
  GLint proji  = glGetUniformLocation(data.meshgl.shader_mesh,"proj");
  glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
  glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
  glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());

  // Light parameters
  GLint specular_exponenti    = glGetUniformLocation(data.meshgl.shader_mesh,"specular_exponent");
  GLint light_position_worldi = glGetUniformLocation(data.meshgl.shader_mesh,"light_position_world");
  GLint lighting_factori      = glGetUniformLocation(data.meshgl.shader_mesh,"lighting_factor");
  GLint fixed_colori          = glGetUniformLocation(data.meshgl.shader_mesh,"fixed_color");
  GLint texture_factori       = glGetUniformLocation(data.meshgl.shader_mesh,"texture_factor");

  glUniform1f(specular_exponenti, data.shininess);
  Vector3f rev_light = -1.*light_position;
  glUniform3fv(light_position_worldi, 1, rev_light.data());
  glUniform1f(lighting_factori, lighting_factor); // enables lighting
  glUniform4f(fixed_colori, 0.0, 0.0, 0.0, 0.0);

  if (data.V.rows()>0)
  {
    // Render fill
    if (data.show_faces)
    {
      // Texture
      glUniform1f(texture_factori, data.show_texture ? 1.0f : 0.0f);
      data.meshgl.draw_mesh(true);
      glUniform1f(texture_factori, 0.0f);
    }

    // Render wireframe
    if (data.show_lines)
    {
      glLineWidth(data.line_width);
      glUniform4f(fixed_colori, 
        data.line_color[0], 
        data.line_color[1],
        data.line_color[2], 1.0f);
      data.meshgl.draw_mesh(false);
      glUniform4f(fixed_colori, 0.0f, 0.0f, 0.0f, 0.0f);
    }
  }

  if (data.show_overlay)
  {
    if (data.show_overlay_depth)
      glEnable(GL_DEPTH_TEST);
    else
      glDisable(GL_DEPTH_TEST);

    if (data.lines.rows() > 0)
    {
      data.meshgl.bind_overlay_lines();
      modeli = glGetUniformLocation(data.meshgl.shader_overlay_lines,"model");
      viewi  = glGetUniformLocation(data.meshgl.shader_overlay_lines,"view");
      proji  = glGetUniformLocation(data.meshgl.shader_overlay_lines,"proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      // This must be enabled, otherwise glLineWidth has no effect
      glEnable(GL_LINE_SMOOTH);
      glLineWidth(data.line_width);

      data.meshgl.draw_overlay_lines();
    }

    if (data.points.rows() > 0)
    {
      data.meshgl.bind_overlay_points();
      modeli = glGetUniformLocation(data.meshgl.shader_overlay_points,"model");
      viewi  = glGetUniformLocation(data.meshgl.shader_overlay_points,"view");
      proji  = glGetUniformLocation(data.meshgl.shader_overlay_points,"proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      glPointSize(data.point_size);

      data.meshgl.draw_overlay_points();
    }

    glEnable(GL_DEPTH_TEST);
  }

}

IGL_INLINE void igl::opengl::ViewerCore::draw_buffer(ViewerData& data,
  bool update_matrices,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A)
{
  assert(R.rows() == G.rows() && G.rows() == B.rows() && B.rows() == A.rows());
  assert(R.cols() == G.cols() && G.cols() == B.cols() && B.cols() == A.cols());

  unsigned x = R.rows();
  unsigned y = R.cols();

  // Create frame buffer
  GLuint frameBuffer;
  glGenFramebuffers(1, &frameBuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

  // Create texture to hold color buffer
  GLuint texColorBuffer;
  glGenTextures(1, &texColorBuffer);
  glBindTexture(GL_TEXTURE_2D, texColorBuffer);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, x, y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);

  // Create Renderbuffer Object to hold depth and stencil buffers
  GLuint rboDepthStencil;
  glGenRenderbuffers(1, &rboDepthStencil);
  glBindRenderbuffer(GL_RENDERBUFFER, rboDepthStencil);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, x, y);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rboDepthStencil);

  assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);

  glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

  // Clear the buffer
  glClearColor(background_color(0), background_color(1), background_color(2), 0.f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Save old viewport
  Eigen::Vector4f viewport_ori = viewport;
  viewport << 0,0,x,y;

  // Draw
  draw(data,update_matrices);

  // Restore viewport
  viewport = viewport_ori;

  // Copy back in the given Eigen matrices
  GLubyte* pixels = (GLubyte*)calloc(x*y*4,sizeof(GLubyte));
  glReadPixels
  (
   0, 0,
   x, y,
   GL_RGBA, GL_UNSIGNED_BYTE, pixels
   );

  int count = 0;
  for (unsigned j=0; j<y; ++j)
  {
    for (unsigned i=0; i<x; ++i)
    {
      R(i,j) = pixels[count*4+0];
      G(i,j) = pixels[count*4+1];
      B(i,j) = pixels[count*4+2];
      A(i,j) = pixels[count*4+3];
      ++count;
    }
  }

  // Clean up
  free(pixels);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glDeleteRenderbuffers(1, &rboDepthStencil);
  glDeleteTextures(1, &texColorBuffer);
  glDeleteFramebuffers(1, &frameBuffer);
}

IGL_INLINE void igl::opengl::ViewerCore::set_rotation_type(
  const igl::opengl::ViewerCore::RotationType & value)
{
  using namespace Eigen;
  using namespace std;
  const RotationType old_rotation_type = rotation_type;
  rotation_type = value;
  if(rotation_type == ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP &&
    old_rotation_type != ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
  {
    snap_to_fixed_up(Quaternionf(trackball_angle),trackball_angle);
  }
}


IGL_INLINE igl::opengl::ViewerCore::ViewerCore()
{
  // Default colors
  background_color << 0.3f, 0.3f, 0.5f, 1.0f;

  // Default lights settings
  light_position << 0.0f, -0.30f, -5.0f;
  lighting_factor = 1.0f; //on

  // Default trackball
  trackball_angle = Eigen::Quaternionf::Identity();
  set_rotation_type(ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);

  // Defalut model viewing parameters
  model_zoom = 1.0f;
  model_translation << 0,0,0;

  // Camera parameters
  camera_zoom = 1.0f;
  orthographic = false;
  camera_view_angle = 45.0;
  camera_dnear = 1.0;
  camera_dfar = 100.0;
  camera_eye << 0, 0, 5;
  camera_center << 0, 0, 0;
  camera_up << 0, 1, 0;

  depth_test = true;

  is_animating = false;
  animation_max_fps = 30.;

  viewport.setZero();
}

IGL_INLINE void igl::opengl::ViewerCore::init()
{
}

IGL_INLINE void igl::opengl::ViewerCore::shut()
{
}
