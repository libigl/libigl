// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "ViewerCore.h"
#include <igl/quat_to_mat.h>
#include <igl/massmatrix.h>
#include <Eigen/Geometry>


Eigen::Matrix4f lookAt (
                        const Eigen::Vector3f& eye,
                        const Eigen::Vector3f& center,
                        const Eigen::Vector3f& up)
{
  Eigen::Vector3f f = (center - eye).normalized();
  Eigen::Vector3f s = f.cross(up).normalized();
  Eigen::Vector3f u = s.cross(f);

  Eigen::Matrix4f Result = Eigen::Matrix4f::Identity();
  Result(0,0) = s(0);
  Result(0,1) = s(1);
  Result(0,2) = s(2);
  Result(1,0) = u(0);
  Result(1,1) = u(1);
  Result(1,2) = u(2);
  Result(2,0) =-f(0);
  Result(2,1) =-f(1);
  Result(2,2) =-f(2);
  Result(0,3) =-s.transpose() * eye;
  Result(1,3) =-u.transpose() * eye;
  Result(2,3) = f.transpose() * eye;
  return Result;
}

Eigen::Matrix4f ortho (
                       const float left,
                       const float right,
                       const float bottom,
                       const float top,
                       const float zNear,
                       const float zFar
                       )
{
  Eigen::Matrix4f Result = Eigen::Matrix4f::Identity();
  Result(0,0) = 2.0f / (right - left);
  Result(1,1) = 2.0f / (top - bottom);
  Result(2,2) = - 2.0f / (zFar - zNear);
  Result(0,3) = - (right + left) / (right - left);
  Result(1,3) = - (top + bottom) / (top - bottom);
  Result(2,3) = - (zFar + zNear) / (zFar - zNear);
  return Result;
}

Eigen::Matrix4f frustum (
                         const float left,
                         const float right,
                         const float bottom,
                         const float top,
                         const float nearVal,
                         const float farVal)
{
  Eigen::Matrix4f Result = Eigen::Matrix4f::Zero();
  Result(0,0) = (2.0f * nearVal) / (right - left);
  Result(1,1) = (2.0f * nearVal) / (top - bottom);
  Result(0,2) = (right + left) / (right - left);
  Result(1,2) = (top + bottom) / (top - bottom);
  Result(2,2) = -(farVal + nearVal) / (farVal - nearVal);
  Result(3,2) = -1.0f;
  Result(2,3) = -(2.0f * farVal * nearVal) / (farVal - nearVal);
  return Result;
}

Eigen::Matrix4f scale (const Eigen::Matrix4f& m,
                       const Eigen::Vector3f& v)
{
  Eigen::Matrix4f Result;
  Result.col(0) = m.col(0).array() * v(0);
  Result.col(1) = m.col(1).array() * v(1);
  Result.col(2) = m.col(2).array() * v(2);
  Result.col(3) = m.col(3);
  return Result;
}

Eigen::Matrix4f translate(
                          const Eigen::Matrix4f& m,
                          const Eigen::Vector3f& v)
{
  Eigen::Matrix4f Result = m;
  Result.col(3) = m.col(0).array() * v(0) + m.col(1).array() * v(1) + m.col(2).array() * v(2) + m.col(3).array();
  return Result;
}



void igl::ViewerCore::InitSerialization()
{
  #ifdef ENABLE_XML_SERIALIZATION
  xmlSerializer->Add(shininess, "shininess");
  xmlSerializer->Add(background_color, "background_color");
  xmlSerializer->Add(line_color, "line_color");
  xmlSerializer->Add(light_position, "light_position");
  xmlSerializer->Add(lighting_factor, "lighting_factor");
  xmlSerializer->Add(trackball_angle, "trackball_angle");
  xmlSerializer->Add(model_zoom, "model_zoom");
  xmlSerializer->Add(model_translation, "model_translation");
  xmlSerializer->Add(model_zoom_uv, "model_zoom_uv");
  xmlSerializer->Add(model_translation_uv, "model_translation_uv");
  xmlSerializer->Add(camera_zoom, "camera_zoom");
  xmlSerializer->Add(orthographic, "orthographic");
  xmlSerializer->Add(camera_eye, "camera_eye");
  xmlSerializer->Add(camera_up, "camera_up");
  xmlSerializer->Add(camera_center, "camera_center");
  xmlSerializer->Add(camera_view_angle, "camera_view_angle");
  xmlSerializer->Add(camera_dnear, "camera_dnear");
  xmlSerializer->Add(camera_dfar, "camera_dfar");
  xmlSerializer->Add(show_overlay, "show_overlay");
  xmlSerializer->Add(show_overlay_depth, "show_overlay_depth");
  xmlSerializer->Add(show_texture, "show_texture");
  xmlSerializer->Add(show_faces, "show_faces");
  xmlSerializer->Add(show_lines, "show_lines");
  xmlSerializer->Add(show_vertid, "show_vertid");
  xmlSerializer->Add(show_faceid, "show_faceid");
  xmlSerializer->Add(point_size, "point_size");
  xmlSerializer->Add(line_width, "line_width");
  xmlSerializer->Add(invert_normals, "invert_normals");
  xmlSerializer->Add(face_based, "face_based");
  xmlSerializer->Add(face_based, "object_scale");
  xmlSerializer->Add(viewport, "viewport");
  xmlSerializer->Add(view, "view");
  xmlSerializer->Add(model, "model");
  xmlSerializer->Add(proj, "proj");

  #endif
}

IGL_INLINE void igl::ViewerCore::align_camera_center(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F)
{
  get_scale_and_shift_to_fit_mesh(V,F,model_zoom,model_translation);
  object_scale = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
}

IGL_INLINE void igl::ViewerCore::get_scale_and_shift_to_fit_mesh(
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
  Eigen::RowVector3d min_point = V.colwise().minCoeff();
  Eigen::RowVector3d max_point = V.colwise().maxCoeff();
  Eigen::RowVector3d centroid  = 0.5*(min_point + max_point);

  shift = -centroid.cast<float>();
  double x_scale = fabs(max_point[0] - min_point[0]);
  double y_scale = fabs(max_point[1] - min_point[1]);
  double z_scale = fabs(max_point[2] - min_point[2]);
  zoom = 2.0 / std::max(z_scale,std::max(x_scale,y_scale));
}

IGL_INLINE void igl::ViewerCore::clear_framebuffers()
{
  glClearColor(background_color[0],
               background_color[1],
               background_color[2],
               1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

IGL_INLINE void igl::ViewerCore::draw(ViewerData& data, OpenGL_state& opengl)
{
  using namespace std;
  using namespace Eigen;

  glEnable(GL_DEPTH_TEST);

  /* Bind and potentially refresh mesh/line/point data */
  if (data.dirty)
  {
    opengl.set_data(data, invert_normals);
    data.dirty = ViewerData::DIRTY_NONE;
  }
  opengl.bind_mesh();

  // Initialize uniform
  glViewport(viewport(0), viewport(1), viewport(2), viewport(3));

  model = Eigen::Matrix4f::Identity();
  view  = Eigen::Matrix4f::Identity();
  proj  = Eigen::Matrix4f::Identity();

  // Set view
  view = lookAt(Eigen::Vector3f(camera_eye[0], camera_eye[1], camera_eye[2]),
                Eigen::Vector3f(camera_center[0], camera_center[1], camera_center[2]),
                Eigen::Vector3f(camera_up[0], camera_up[1], camera_up[2]));

  float width  = viewport(2);
  float height = viewport(3);

  // Set projection
  if (orthographic)
  {
    float length = (camera_eye - camera_center).norm();
    float h = tan(camera_view_angle/360.0 * M_PI) * (length);
    proj = ortho(-h*width/height, h*width/height, -h, h, camera_dnear, camera_dfar);
  }
  else
  {
    float fH = tan(camera_view_angle / 360.0 * M_PI) * camera_dnear;
    float fW = fH * (double)width/(double)height;
    proj = frustum(-fW, fW, -fH, fH, camera_dnear, camera_dfar);
  }
  // end projection

  // Set model transformation
  float mat[16];
  igl::quat_to_mat(trackball_angle.data(), mat);

  for (unsigned i=0;i<4;++i)
    for (unsigned j=0;j<4;++j)
      model(i,j) = mat[i+4*j];

  model = scale(model, Eigen::Vector3f(camera_zoom,camera_zoom,camera_zoom));
  model = scale(model, Eigen::Vector3f(model_zoom,model_zoom,model_zoom));
  model = translate(model, Eigen::Vector3f(model_translation[0],model_translation[1],model_translation[2]));

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
        textrenderer.DrawText(data.V.row(i), data.V_normals.row(i), to_string(i));
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

IGL_INLINE igl::ViewerCore::ViewerCore()
#ifdef ENABLE_XML_SERIALIZATION
: XMLSerialization("Core")
#endif
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
  trackball_angle << 0.0f, 0.0f, 0.0f, 1.0f;

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

  // Default point size / line width
  point_size = 15;
  line_width = 0.5f;
  is_animating = false;
  animation_max_fps = 30.;
}

IGL_INLINE void igl::ViewerCore::init()
{
  textrenderer.Init();
}

IGL_INLINE void igl::ViewerCore::shut()
{
  textrenderer.Shut();
}
