// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "OpenGL_state.h"
#include "ViewerData.h"

IGL_INLINE void igl::viewer::OpenGL_state::init_buffers()
{
  // Mesh: Vertex Array Object & Buffer objects
  glGenVertexArrays(1, &vao_mesh);
  glBindVertexArray(vao_mesh);
  glGenBuffers(1, &vbo_V);
  glGenBuffers(1, &vbo_V_normals);
  glGenBuffers(1, &vbo_V_ambient);
  glGenBuffers(1, &vbo_V_diffuse);
  glGenBuffers(1, &vbo_V_specular);
  glGenBuffers(1, &vbo_V_uv);
  glGenBuffers(1, &vbo_F);
  glGenTextures(1, &vbo_tex);

  // Line overlay
  glGenVertexArrays(1, &vao_overlay_lines);
  glBindVertexArray(vao_overlay_lines);
  glGenBuffers(1, &vbo_lines_F);
  glGenBuffers(1, &vbo_lines_V);
  glGenBuffers(1, &vbo_lines_V_colors);

  // Point overlay
  glGenVertexArrays(1, &vao_overlay_points);
  glBindVertexArray(vao_overlay_points);
  glGenBuffers(1, &vbo_points_F);
  glGenBuffers(1, &vbo_points_V);
  glGenBuffers(1, &vbo_points_V_colors);

  dirty = ViewerData::DIRTY_ALL;
}

IGL_INLINE void igl::viewer::OpenGL_state::free_buffers()
{
  glDeleteVertexArrays(1, &vao_mesh);
  glDeleteVertexArrays(1, &vao_overlay_lines);
  glDeleteVertexArrays(1, &vao_overlay_points);

  glDeleteBuffers(1, &vbo_V);
  glDeleteBuffers(1, &vbo_V_normals);
  glDeleteBuffers(1, &vbo_V_ambient);
  glDeleteBuffers(1, &vbo_V_diffuse);
  glDeleteBuffers(1, &vbo_V_specular);
  glDeleteBuffers(1, &vbo_V_uv);
  glDeleteBuffers(1, &vbo_F);
  glDeleteBuffers(1, &vbo_lines_F);
  glDeleteBuffers(1, &vbo_lines_V);
  glDeleteBuffers(1, &vbo_lines_V_colors);
  glDeleteBuffers(1, &vbo_points_F);
  glDeleteBuffers(1, &vbo_points_V);
  glDeleteBuffers(1, &vbo_points_V_colors);

  glDeleteTextures(1, &vbo_tex);
}

IGL_INLINE void igl::viewer::OpenGL_state::set_data(const igl::viewer::ViewerData &data, bool invert_normals)
{
  bool per_corner_uv = (data.F_uv.rows() == data.F.rows());
  bool per_corner_normals = (data.F_normals.rows() == 3 * data.F.rows());

  dirty |= data.dirty;

  if (!data.face_based)
  {
    if (!per_corner_uv)
    {
      // Vertex positions
      if (dirty & ViewerData::DIRTY_POSITION)
        V_vbo = (data.V.transpose()).cast<float>();

      // Vertex normals
      if (dirty & ViewerData::DIRTY_NORMAL)
      {
        V_normals_vbo = (data.V_normals.transpose()).cast<float>();
        if (invert_normals)
          V_normals_vbo = -V_normals_vbo;
      }

      // Per-vertex material settings
      if (dirty & ViewerData::DIRTY_AMBIENT)
        V_ambient_vbo = (data.V_material_ambient.transpose()).cast<float>();
      if (dirty & ViewerData::DIRTY_DIFFUSE)
        V_diffuse_vbo = (data.V_material_diffuse.transpose()).cast<float>();
      if (dirty & ViewerData::DIRTY_SPECULAR)
        V_specular_vbo = (data.V_material_specular.transpose()).cast<float>();

      // Face indices
      if (dirty & ViewerData::DIRTY_FACE)
        F_vbo = (data.F.transpose()).cast<unsigned>();

      // Texture coordinates
      if (dirty & ViewerData::DIRTY_UV)
        V_uv_vbo = (data.V_uv.transpose()).cast<float>();
    }
    else
    {
      // Per vertex properties with per corner UVs
      if (dirty & ViewerData::DIRTY_POSITION)
      {
        V_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_vbo.col(i*3+j) = data.V.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & ViewerData::DIRTY_AMBIENT)
      {
        V_ambient_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_ambient_vbo.col (i*3+j) = data.V_material_ambient.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & ViewerData::DIRTY_DIFFUSE)
      {
        V_diffuse_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_diffuse_vbo.col (i*3+j) = data.V_material_diffuse.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & ViewerData::DIRTY_SPECULAR)
      {
        V_specular_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_specular_vbo.col(i*3+j) = data.V_material_specular.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & ViewerData::DIRTY_NORMAL)
      {
        V_normals_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_normals_vbo.col (i*3+j) = data.V_normals.row(data.F(i,j)).transpose().cast<float>();

        if (invert_normals)
          V_normals_vbo = -V_normals_vbo;
      }

      if (dirty & ViewerData::DIRTY_FACE)
      {
        F_vbo.resize(3,data.F.rows());
        for (unsigned i=0; i<data.F.rows();++i)
          F_vbo.col(i) << i*3+0, i*3+1, i*3+2;
      }

      if (dirty & ViewerData::DIRTY_UV)
      {
        V_uv_vbo.resize(2,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_uv_vbo.col(i*3+j) = data.V_uv.row(data.F(i,j)).transpose().cast<float>();
      }
    }
  }
  else
  {
    if (dirty & ViewerData::DIRTY_POSITION)
    {
      V_vbo.resize(3,data.F.rows()*3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_vbo.col(i*3+j) = data.V.row(data.F(i,j)).transpose().cast<float>();
    }

    if (dirty & ViewerData::DIRTY_AMBIENT)
    {
      V_ambient_vbo.resize(3,data.F.rows()*3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_ambient_vbo.col (i*3+j) = data.F_material_ambient.row(i).transpose().cast<float>();
    }

    if (dirty & ViewerData::DIRTY_DIFFUSE)
    {
      V_diffuse_vbo.resize(3,data.F.rows()*3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_diffuse_vbo.col (i*3+j) = data.F_material_diffuse.row(i).transpose().cast<float>();
    }

    if (dirty & ViewerData::DIRTY_SPECULAR)
    {
      V_specular_vbo.resize(3,data.F.rows()*3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_specular_vbo.col(i*3+j) = data.F_material_specular.row(i).transpose().cast<float>();
    }

    if (dirty & ViewerData::DIRTY_NORMAL)
    {
      V_normals_vbo.resize(3,data.F.rows()*3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_normals_vbo.col (i*3+j) =
             per_corner_normals ?
               data.F_normals.row(i*3+j).transpose().cast<float>() :
               data.F_normals.row(i).transpose().cast<float>();

      if (invert_normals)
        V_normals_vbo = -V_normals_vbo;
    }

    if (dirty & ViewerData::DIRTY_FACE)
    {
      F_vbo.resize(3,data.F.rows());
      for (unsigned i=0; i<data.F.rows();++i)
        F_vbo.col(i) << i*3+0, i*3+1, i*3+2;
    }

    if (dirty & ViewerData::DIRTY_UV)
    {
        V_uv_vbo.resize(2,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_uv_vbo.col(i*3+j) = data.V_uv.row(per_corner_uv ? data.F_uv(i,j) : data.F(i,j)).transpose().cast<float>();
    }
  }

  if (dirty & ViewerData::DIRTY_TEXTURE)
  {
    tex_u = data.texture_R.rows();
    tex_v = data.texture_R.cols();
    tex.resize(data.texture_R.size()*3);
    for (unsigned i=0;i<data.texture_R.size();++i)
    {
      tex(i*3+0) = data.texture_R(i);
      tex(i*3+1) = data.texture_G(i);
      tex(i*3+2) = data.texture_B(i);
    }
  }

  if (dirty & ViewerData::DIRTY_OVERLAY_LINES)
  {
    lines_V_vbo.resize(3, data.lines.rows()*2);
    lines_V_colors_vbo.resize(3, data.lines.rows()*2);
    lines_F_vbo.resize(1, data.lines.rows()*2);
    for (unsigned i=0; i<data.lines.rows();++i)
    {
      lines_V_vbo.col(2*i+0) = data.lines.block<1, 3>(i, 0).transpose().cast<float>();
      lines_V_vbo.col(2*i+1) = data.lines.block<1, 3>(i, 3).transpose().cast<float>();
      lines_V_colors_vbo.col(2*i+0) = data.lines.block<1, 3>(i, 6).transpose().cast<float>();
      lines_V_colors_vbo.col(2*i+1) = data.lines.block<1, 3>(i, 6).transpose().cast<float>();
      lines_F_vbo(2*i+0) = 2*i+0;
      lines_F_vbo(2*i+1) = 2*i+1;
    }
  }

  if (dirty & ViewerData::DIRTY_OVERLAY_POINTS)
  {
    points_V_vbo.resize(3, data.points.rows());
    points_V_colors_vbo.resize(3, data.points.rows());
    points_F_vbo.resize(1, data.points.rows());
    for (unsigned i=0; i<data.points.rows();++i)
    {
      points_V_vbo.col(i) = data.points.block<1, 3>(i, 0).transpose().cast<float>();
      points_V_colors_vbo.col(i) = data.points.block<1, 3>(i, 3).transpose().cast<float>();
      points_F_vbo(i) = i;
    }
  }
}

IGL_INLINE void igl::viewer::OpenGL_state::bind_mesh()
{
  glBindVertexArray(vao_mesh);
  shader_mesh.bind();
  shader_mesh.bindVertexAttribArray("position", vbo_V, V_vbo, dirty & ViewerData::DIRTY_POSITION);
  shader_mesh.bindVertexAttribArray("normal", vbo_V_normals, V_normals_vbo, dirty & ViewerData::DIRTY_NORMAL);
  shader_mesh.bindVertexAttribArray("Ka", vbo_V_ambient, V_ambient_vbo, dirty & ViewerData::DIRTY_AMBIENT);
  shader_mesh.bindVertexAttribArray("Kd", vbo_V_diffuse, V_diffuse_vbo, dirty & ViewerData::DIRTY_DIFFUSE);
  shader_mesh.bindVertexAttribArray("Ks", vbo_V_specular, V_specular_vbo, dirty & ViewerData::DIRTY_SPECULAR);
  shader_mesh.bindVertexAttribArray("texcoord", vbo_V_uv, V_uv_vbo, dirty & ViewerData::DIRTY_UV);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_F);
  if (dirty & ViewerData::DIRTY_FACE)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*F_vbo.size(), F_vbo.data(), GL_DYNAMIC_DRAW);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, vbo_tex);
  if (dirty & ViewerData::DIRTY_TEXTURE)
  {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_u, tex_v, 0, GL_RGB, GL_UNSIGNED_BYTE, tex.data());
  }
  glUniform1i(shader_mesh.uniform("tex"), 0);
  dirty &= ~ViewerData::DIRTY_MESH;
}

IGL_INLINE void igl::viewer::OpenGL_state::bind_overlay_lines()
{
  bool is_dirty = dirty & ViewerData::DIRTY_OVERLAY_LINES;

  glBindVertexArray(vao_overlay_lines);
  shader_overlay_lines.bind();
  shader_overlay_lines.bindVertexAttribArray("position", vbo_lines_V, lines_V_vbo, is_dirty);
  shader_overlay_lines.bindVertexAttribArray("color", vbo_lines_V_colors, lines_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_lines_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*lines_F_vbo.size(), lines_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~ViewerData::DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::viewer::OpenGL_state::bind_overlay_points()
{
  bool is_dirty = dirty & ViewerData::DIRTY_OVERLAY_POINTS;

  glBindVertexArray(vao_overlay_points);
  shader_overlay_points.bind();
  shader_overlay_points.bindVertexAttribArray("position", vbo_points_V, points_V_vbo, is_dirty);
  shader_overlay_points.bindVertexAttribArray("color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*points_F_vbo.size(), points_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~ViewerData::DIRTY_OVERLAY_POINTS;
}

IGL_INLINE void igl::viewer::OpenGL_state::draw_mesh(bool solid)
{
  glPolygonMode(GL_FRONT_AND_BACK, solid ? GL_FILL : GL_LINE);

  /* Avoid Z-buffer fighting between filled triangles & wireframe lines */
  if (solid)
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
  }
  glDrawElements(GL_TRIANGLES, 3*F_vbo.cols(), GL_UNSIGNED_INT, 0);

  glDisable(GL_POLYGON_OFFSET_FILL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

IGL_INLINE void igl::viewer::OpenGL_state::draw_overlay_lines()
{
  glDrawElements(GL_LINES, lines_F_vbo.cols(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::viewer::OpenGL_state::draw_overlay_points()
{
  glDrawElements(GL_POINTS, points_F_vbo.cols(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::viewer::OpenGL_state::init()
{
  std::string mesh_vertex_shader_string =
  "#version 150\n"
  "uniform mat4 model;"
  "uniform mat4 view;"
  "uniform mat4 proj;"
  "in vec3 position;"
  "in vec3 normal;"
  "out vec3 position_eye;"
  "out vec3 normal_eye;"
  "in vec3 Ka;"
  "in vec3 Kd;"
  "in vec3 Ks;"
  "in vec2 texcoord;"
  "out vec2 texcoordi;"
  "out vec3 Kai;"
  "out vec3 Kdi;"
  "out vec3 Ksi;"

  "void main()"
  "{"
  "  position_eye = vec3 (view * model * vec4 (position, 1.0));"
  "  normal_eye = vec3 (view * model * vec4 (normal, 0.0));"
  "  normal_eye = normalize(normal_eye);"
  "  gl_Position = proj * vec4 (position_eye, 1.0);" //proj * view * model * vec4(position, 1.0);"
  "  Kai = Ka;"
  "  Kdi = Kd;"
  "  Ksi = Ks;"
  "  texcoordi = texcoord;"
  "}";

  std::string mesh_fragment_shader_string =
  "#version 150\n"
  "uniform mat4 model;"
  "uniform mat4 view;"
  "uniform mat4 proj;"
  "uniform vec4 fixed_color;"
  "in vec3 position_eye;"
  "in vec3 normal_eye;"
  "uniform vec3 light_position_world;"
  "vec3 Ls = vec3 (1, 1, 1);"
  "vec3 Ld = vec3 (1, 1, 1);"
  "vec3 La = vec3 (1, 1, 1);"
  "in vec3 Ksi;"
  "in vec3 Kdi;"
  "in vec3 Kai;"
  "in vec2 texcoordi;"
  "uniform sampler2D tex;"
  "uniform float specular_exponent;"
  "uniform float lighting_factor;"
  "uniform float texture_factor;"
  "out vec4 outColor;"
  "void main()"
  "{"
  "vec3 Ia = La * Kai;"    // ambient intensity

  "vec3 light_position_eye = vec3 (view * vec4 (light_position_world, 1.0));"
  "vec3 vector_to_light_eye = light_position_eye - position_eye;"
  "vec3 direction_to_light_eye = normalize (vector_to_light_eye);"
  "float dot_prod = dot (direction_to_light_eye, normal_eye);"
  "float clamped_dot_prod = max (dot_prod, 0.0);"
  "vec3 Id = Ld * Kdi * clamped_dot_prod;"    // Diffuse intensity

  "vec3 reflection_eye = reflect (-direction_to_light_eye, normal_eye);"
  "vec3 surface_to_viewer_eye = normalize (-position_eye);"
  "float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);"
  "dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);"
  "float specular_factor = pow (dot_prod_specular, specular_exponent);"
  "vec3 Is = Ls * Ksi * specular_factor;"    // specular intensity
  "vec4 color = vec4(lighting_factor * (Is + Id) + Ia, 1.0) + vec4((1.0-lighting_factor) * Kdi,1.0);"
  "outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;"
  "if (fixed_color != vec4(0.0)) outColor = fixed_color;"
  "}";

  std::string overlay_vertex_shader_string =
  "#version 150\n"
  "uniform mat4 model;"
  "uniform mat4 view;"
  "uniform mat4 proj;"
  "in vec3 position;"
  "in vec3 color;"
  "out vec3 color_frag;"

  "void main()"
  "{"
  "  gl_Position = proj * view * model * vec4 (position, 1.0);"
  "  color_frag = color;"
  "}";

  std::string overlay_fragment_shader_string =
  "#version 150\n"
  "in vec3 color_frag;"
  "out vec4 outColor;"
  "void main()"
  "{"
  "  outColor = vec4(color_frag, 1.0);"
  "}";

  std::string overlay_point_fragment_shader_string =
  "#version 150\n"
  "in vec3 color_frag;"
  "out vec4 outColor;"
  "void main()"
  "{"
  "  if (length(gl_PointCoord - vec2(0.5)) > 0.5)"
  "    discard;"
  "  outColor = vec4(color_frag, 1.0);"
  "}";

  init_buffers();
  shader_mesh.init(mesh_vertex_shader_string,
      mesh_fragment_shader_string, "outColor");
  shader_overlay_lines.init(overlay_vertex_shader_string,
      overlay_fragment_shader_string, "outColor");
  shader_overlay_points.init(overlay_vertex_shader_string,
      overlay_point_fragment_shader_string, "outColor");
}

IGL_INLINE void igl::viewer::OpenGL_state::free()
{
  shader_mesh.free();
  shader_overlay_lines.free();
  shader_overlay_points.free();
  free_buffers();
}
