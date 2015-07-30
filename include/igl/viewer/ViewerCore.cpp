// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerCore.h"
#include <igl/quat_to_mat.h>
#include <igl/look_at.h>
#include <igl/frustum.h>
#include <igl/ortho.h>
#include <igl/massmatrix.h>
#include <igl/barycenter.h>
#include <Eigen/Geometry>
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <igl/serialize.h>
namespace igl {
  namespace serialization {

    IGL_INLINE void serialization(bool s,igl::viewer::ViewerCore& obj,std::vector<char>& buffer)
    {
      SERIALIZE_MEMBER(shininess);

      SERIALIZE_MEMBER(background_color);
      SERIALIZE_MEMBER(line_color);

      SERIALIZE_MEMBER(light_position);
      SERIALIZE_MEMBER(lighting_factor);

      SERIALIZE_MEMBER(trackball_angle);

      SERIALIZE_MEMBER(model_zoom);
      SERIALIZE_MEMBER(model_translation);

      SERIALIZE_MEMBER(model_zoom_uv);
      SERIALIZE_MEMBER(model_translation_uv);

      SERIALIZE_MEMBER(object_scale);

      SERIALIZE_MEMBER(camera_zoom);
      SERIALIZE_MEMBER(orthographic);
      SERIALIZE_MEMBER(camera_view_angle);
      SERIALIZE_MEMBER(camera_dnear);
      SERIALIZE_MEMBER(camera_dfar);
      SERIALIZE_MEMBER(camera_eye);
      SERIALIZE_MEMBER(camera_center);
      SERIALIZE_MEMBER(camera_up);

      SERIALIZE_MEMBER(show_faces);
      SERIALIZE_MEMBER(show_lines);
      SERIALIZE_MEMBER(invert_normals);
      SERIALIZE_MEMBER(show_overlay);
      SERIALIZE_MEMBER(show_overlay_depth);
      SERIALIZE_MEMBER(show_vertid);
      SERIALIZE_MEMBER(show_faceid);
      SERIALIZE_MEMBER(show_texture);
      SERIALIZE_MEMBER(depth_test);

      SERIALIZE_MEMBER(point_size);
      SERIALIZE_MEMBER(line_width);
      SERIALIZE_MEMBER(is_animating);
      SERIALIZE_MEMBER(animation_max_fps);

      SERIALIZE_MEMBER(viewport);
      SERIALIZE_MEMBER(view);
      SERIALIZE_MEMBER(model);
      SERIALIZE_MEMBER(proj);
    }

    IGL_INLINE void serialize(const igl::viewer::ViewerCore& obj,std::vector<char>& buffer)
    {
      serialization(true,const_cast<igl::viewer::ViewerCore&>(obj),buffer);
    }

    IGL_INLINE void deserialize(igl::viewer::ViewerCore& obj,const std::vector<char>& buffer)
    {
      serialization(false,obj,const_cast<std::vector<char>&>(buffer));
    }
  }
}
#endif

IGL_INLINE void igl::viewer::ViewerCore::align_camera_center(
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

IGL_INLINE void igl::viewer::ViewerCore::get_scale_and_shift_to_fit_mesh(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  float& zoom,
  Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  //Eigen::SparseMatrix<double> M;
  //igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  //const auto & MV = M*V;
  //Eigen::RowVector3d centroid  = MV.colwise().sum()/M.diagonal().sum();
  Eigen::MatrixXd BC;
  igl::barycenter(V,F,BC);
  Eigen::RowVector3d min_point = BC.colwise().minCoeff();
  Eigen::RowVector3d max_point = BC.colwise().maxCoeff();
  Eigen::RowVector3d centroid  = 0.5*(min_point + max_point);

  shift = -centroid.cast<float>();
  double x_scale = fabs(max_point[0] - min_point[0]);
  double y_scale = fabs(max_point[1] - min_point[1]);
  double z_scale = fabs(max_point[2] - min_point[2]);
  zoom = 2.0 / std::max(z_scale,std::max(x_scale,y_scale));
}

IGL_INLINE void igl::viewer::ViewerCore::align_camera_center(
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

IGL_INLINE void igl::viewer::ViewerCore::get_scale_and_shift_to_fit_mesh(
                                                                 const Eigen::MatrixXd& V,
                                                                 float& zoom,
                                                                 Eigen::Vector3f& shift)
{
  if (V.rows() == 0)
    return;

  //Eigen::SparseMatrix<double> M;
  //igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  //const auto & MV = M*V;
  //Eigen::RowVector3d centroid  = MV.colwise().sum()/M.diagonal().sum();
  Eigen::RowVector3d min_point = V.colwise().minCoeff();
  Eigen::RowVector3d max_point = V.colwise().maxCoeff();
  Eigen::RowVector3d centroid  = 0.5*(min_point + max_point);

  shift = -centroid.cast<float>();
  double x_scale = fabs(max_point[0] - min_point[0]);
  double y_scale = fabs(max_point[1] - min_point[1]);
  double z_scale = fabs(max_point[2] - min_point[2]);
  zoom = 2.0 / std::max(z_scale,std::max(x_scale,y_scale));
}


IGL_INLINE void igl::viewer::ViewerCore::clear_framebuffers()
{
  glClearColor(background_color[0],
               background_color[1],
               background_color[2],
               1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

IGL_INLINE void igl::viewer::ViewerCore::draw(ViewerData& data, OpenGL_state& opengl, bool update_matrices)
{
  using namespace std;
  using namespace Eigen;

  if (depth_test)
    glEnable(GL_DEPTH_TEST);
  else
    glDisable(GL_DEPTH_TEST);

  /* Bind and potentially refresh mesh/line/point data */
  if (data.dirty)
  {
    opengl.set_data(data, invert_normals);
    data.dirty = ViewerData::DIRTY_NONE;
  }
  opengl.bind_mesh();

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
      float h = tan(camera_view_angle/360.0 * M_PI) * (length);
      ortho(-h*width/height, h*width/height, -h, h, camera_dnear, camera_dfar,proj);
    }
    else
    {
      float fH = tan(camera_view_angle / 360.0 * M_PI) * camera_dnear;
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
  GLint modeli = opengl.shader_mesh.uniform("model");
  GLint viewi  = opengl.shader_mesh.uniform("view");
  GLint proji  = opengl.shader_mesh.uniform("proj");
  glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
  glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
  glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());

  // Light parameters
  GLint specular_exponenti    = opengl.shader_mesh.uniform("specular_exponent");
  GLint light_position_worldi = opengl.shader_mesh.uniform("light_position_world");
  GLint lighting_factori      = opengl.shader_mesh.uniform("lighting_factor");
  GLint fixed_colori          = opengl.shader_mesh.uniform("fixed_color");
  GLint texture_factori       = opengl.shader_mesh.uniform("texture_factor");

  glUniform1f(specular_exponenti, shininess);
  Vector3f rev_light = -1.*light_position;
  glUniform3fv(light_position_worldi, 1, rev_light.data());
  glUniform1f(lighting_factori, lighting_factor); // enables lighting
  glUniform4f(fixed_colori, 0.0, 0.0, 0.0, 0.0);

  if (data.V.rows()>0)
  {
    // Render fill
    if (show_faces)
    {
      // Texture
      glUniform1f(texture_factori, show_texture ? 1.0f : 0.0f);
      opengl.draw_mesh(true);
      glUniform1f(texture_factori, 0.0f);
    }

    // Render wireframe
    if (show_lines)
    {
      glLineWidth(line_width);
      glUniform4f(fixed_colori, line_color[0], line_color[1],
        line_color[2], 1.0f);
      opengl.draw_mesh(false);
      glUniform4f(fixed_colori, 0.0f, 0.0f, 0.0f, 0.0f);
    }

    if (show_vertid)
    {
      textrenderer.BeginDraw(view*model, proj, viewport, object_scale);
      for (int i=0; i<data.V.rows(); ++i)
        textrenderer.DrawText(data.V.row(i),data.V_normals.row(i),to_string(i));
      textrenderer.EndDraw();
    }

    if (show_faceid)
    {
      textrenderer.BeginDraw(view*model, proj, viewport, object_scale);

      for (int i=0; i<data.F.rows(); ++i)
      {
        Eigen::RowVector3d p = Eigen::RowVector3d::Zero();
        for (int j=0;j<data.F.cols();++j)
          p += data.V.row(data.F(i,j));
        p /= data.F.cols();

        textrenderer.DrawText(p, data.F_normals.row(i), to_string(i));
      }
      textrenderer.EndDraw();
    }
  }

  if (show_overlay)
  {
    if (show_overlay_depth)
      glEnable(GL_DEPTH_TEST);
    else
      glDisable(GL_DEPTH_TEST);

    if (data.lines.rows() > 0)
    {
      opengl.bind_overlay_lines();
      modeli = opengl.shader_overlay_lines.uniform("model");
      viewi  = opengl.shader_overlay_lines.uniform("view");
      proji  = opengl.shader_overlay_lines.uniform("proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      // This must be enabled, otherwise glLineWidth has no effect
      glEnable(GL_LINE_SMOOTH);
      glLineWidth(line_width);

      opengl.draw_overlay_lines();
    }

    if (data.points.rows() > 0)
    {
      opengl.bind_overlay_points();
      modeli = opengl.shader_overlay_points.uniform("model");
      viewi  = opengl.shader_overlay_points.uniform("view");
      proji  = opengl.shader_overlay_points.uniform("proj");

      glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
      glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
      glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
      glPointSize(point_size);

      opengl.draw_overlay_points();
    }

    if (data.labels_positions.rows() > 0)
    {
      textrenderer.BeginDraw(view*model, proj, viewport, object_scale);
      for (int i=0; i<data.labels_positions.rows(); ++i)
        textrenderer.DrawText(data.labels_positions.row(i), Eigen::Vector3d(0.0,0.0,0.0),
            data.labels_strings[i]);
      textrenderer.EndDraw();
    }

    glEnable(GL_DEPTH_TEST);
  }

}

IGL_INLINE void igl::viewer::ViewerCore::draw_buffer(ViewerData& data,
  OpenGL_state& opengl,
  bool update_matrices,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A)
{
  assert(R.rows() == G.rows() && G.rows() == B.rows() && B.rows() == A.rows());
  assert(R.cols() == G.cols() && G.cols() == B.cols() && B.cols() == A.cols());

  int x = R.rows();
  int y = R.cols();

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
  draw(data,opengl,update_matrices);

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


IGL_INLINE igl::viewer::ViewerCore::ViewerCore()
{
  // Default shininess
  shininess = 35.0f;

  // Default colors
  background_color << 0.3f, 0.3f, 0.5f;
  line_color << 0.0f, 0.0f, 0.0f;

  // Default lights settings
  light_position << 0.0f, -0.30f, -5.0f;
  lighting_factor = 1.0f; //on

  // Default trackball
  trackball_angle = Eigen::Quaternionf::Identity();

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

  // Default visualization options
  show_faces = true;
  show_lines = true;
  invert_normals = false;
  show_overlay = true;
  show_overlay_depth = true;
  show_vertid = false;
  show_faceid = false;
  show_texture = false;
  depth_test = true;

  // Default point size / line width
  point_size = 30;
  line_width = 0.5f;
  is_animating = false;
  animation_max_fps = 30.;
}

IGL_INLINE void igl::viewer::ViewerCore::init()
{
  textrenderer.Init();
}

IGL_INLINE void igl::viewer::ViewerCore::shut()
{
  textrenderer.Shut();
}
