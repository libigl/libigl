// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "MeshGL.h"
#include "bind_vertex_attrib_array.h"
#include "ViewerData.h"
#include "create_shader_program.h"
#include "destroy_shader_program.h"

IGL_INLINE void igl::opengl::MeshGL::init_buffers()
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

IGL_INLINE void igl::opengl::MeshGL::free_buffers()
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

IGL_INLINE void igl::opengl::MeshGL::set_data(
  const igl::opengl::ViewerData &data, 
  bool invert_normals)
{
  bool per_corner_uv = (data.F_uv.rows() == data.F.rows());
  bool per_corner_normals = (data.F_normals.rows() == 3 * data.F.rows());

  dirty |= data.dirty;

  // Input:
  //   X  #F by dim quantity
  // Output:
  //   X_vbo  #F*3 by dim scattering per corner
  const auto per_face = [&data](
      const Eigen::MatrixXd & X,
      RowMatrixXf & X_vbo)
  {
    X_vbo.resize(data.F.rows()*3,3);
    for (unsigned i=0; i<data.F.rows();++i)
      for (unsigned j=0;j<3;++j)
        X_vbo.row(i*3+j) = X.row(i).cast<float>().head(3);
  };

  // Input:
  //   X  #V by dim quantity
  // Output:
  //   X_vbo  #F*3 by dim scattering per corner
  const auto per_corner = [&data](
      const Eigen::MatrixXd & X,
      RowMatrixXf & X_vbo)
  {
    X_vbo.resize(data.F.rows()*3,3);
    for (unsigned i=0; i<data.F.rows();++i)
      for (unsigned j=0;j<3;++j)
        X_vbo.row(i*3+j) = X.row(data.F(i,j)).cast<float>();
  };

  if (!data.face_based)
  {
    if (!(per_corner_uv || per_corner_normals))
    {
      // Vertex positions
      if (dirty & ViewerData::DIRTY_POSITION)
        V_vbo = data.V.cast<float>();

      // Vertex normals
      if (dirty & ViewerData::DIRTY_NORMAL)
      {
        V_normals_vbo = data.V_normals.cast<float>();
        if (invert_normals)
          V_normals_vbo = -V_normals_vbo;
      }

      // Per-vertex material settings
      if (dirty & ViewerData::DIRTY_AMBIENT)
        V_ambient_vbo = data.V_material_ambient.cast<float>();
      if (dirty & ViewerData::DIRTY_DIFFUSE)
        V_diffuse_vbo = data.V_material_diffuse.cast<float>();
      if (dirty & ViewerData::DIRTY_SPECULAR)
        V_specular_vbo = data.V_material_specular.cast<float>();

      // Face indices
      if (dirty & ViewerData::DIRTY_FACE)
        F_vbo = data.F.cast<unsigned>();

      // Texture coordinates
      if (dirty & ViewerData::DIRTY_UV)
        V_uv_vbo = data.V_uv.cast<float>();
    }
    else
    {

      // Per vertex properties with per corner UVs
      if (dirty & ViewerData::DIRTY_POSITION)
      {
        per_corner(data.V,V_vbo);
      }

      if (dirty & ViewerData::DIRTY_AMBIENT)
      {
        V_ambient_vbo.resize(4,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_ambient_vbo.col (i*3+j) = data.V_material_ambient.row(data.F(i,j)).transpose().cast<float>();
      }
      if (dirty & ViewerData::DIRTY_DIFFUSE)
      {
        V_diffuse_vbo.resize(4,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_diffuse_vbo.col (i*3+j) = data.V_material_diffuse.row(data.F(i,j)).transpose().cast<float>();
      }
      if (dirty & ViewerData::DIRTY_SPECULAR)
      {
        V_specular_vbo.resize(4,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_specular_vbo.col(i*3+j) = data.V_material_specular.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & ViewerData::DIRTY_NORMAL)
      {
        V_normals_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
          
            V_normals_vbo.col (i*3+j) = 
                         per_corner_normals ?
               data.F_normals.row(i*3+j).transpose().cast<float>() :
               data.V_normals.row(data.F(i,j)).transpose().cast<float>();


        if (invert_normals)
          V_normals_vbo = -V_normals_vbo;
      }

      if (dirty & ViewerData::DIRTY_FACE)
      {
        F_vbo.resize(data.F.rows(),3);
        for (unsigned i=0; i<data.F.rows();++i)
          F_vbo.row(i) << i*3+0, i*3+1, i*3+2;
      }

      if (dirty & ViewerData::DIRTY_UV)
      {
        V_uv_vbo.resize(data.F.rows()*3,2);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_uv_vbo.row(i*3+j) = 
              data.V_uv.row(per_corner_uv ? 
                data.F_uv(i,j) : data.F(i,j)).cast<float>();
      }
    }
  }
  else
  {
    if (dirty & ViewerData::DIRTY_POSITION)
    {
      per_corner(data.V,V_vbo);
    }

    if (dirty & ViewerData::DIRTY_AMBIENT)
    {
      per_face(data.F_material_ambient,V_ambient_vbo);
    }
    if (dirty & ViewerData::DIRTY_DIFFUSE)
    {
      per_face(data.F_material_diffuse,V_diffuse_vbo);
    }
    if (dirty & ViewerData::DIRTY_SPECULAR)
    {
      per_face(data.F_material_specular,V_specular_vbo);
    }

    if (dirty & ViewerData::DIRTY_NORMAL)
    {
      V_normals_vbo.resize(data.F.rows()*3,3);
      for (unsigned i=0; i<data.F.rows();++i)
        for (unsigned j=0;j<3;++j)
          V_normals_vbo.row(i*3+j) =
             per_corner_normals ?
               data.F_normals.row(i*3+j).cast<float>() :
               data.F_normals.row(i).cast<float>();

      if (invert_normals)
        V_normals_vbo = -V_normals_vbo;
    }

    if (dirty & ViewerData::DIRTY_FACE)
    {
      F_vbo.resize(data.F.rows(),3);
      for (unsigned i=0; i<data.F.rows();++i)
        F_vbo.row(i) << i*3+0, i*3+1, i*3+2;
    }

    if (dirty & ViewerData::DIRTY_UV)
    {
        V_uv_vbo.resize(data.F.rows()*3,2);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_uv_vbo.row(i*3+j) = data.V_uv.row(per_corner_uv ? data.F_uv(i,j) : data.F(i,j)).cast<float>();
    }
  }

  if (dirty & ViewerData::DIRTY_TEXTURE)
  {
    tex_u = data.texture_R.rows();
    tex_v = data.texture_R.cols();
    tex.resize(data.texture_R.size()*4);
    for (unsigned i=0;i<data.texture_R.size();++i)
    {
      tex(i*4+0) = data.texture_R(i);
      tex(i*4+1) = data.texture_G(i);
      tex(i*4+2) = data.texture_B(i);
      tex(i*4+3) = data.texture_A(i);
    }
  }

  if (dirty & ViewerData::DIRTY_OVERLAY_LINES)
  {
    lines_V_vbo.resize(data.lines.rows()*2,3);
    lines_V_colors_vbo.resize(data.lines.rows()*2,3);
    lines_F_vbo.resize(data.lines.rows()*2,1);
    for (unsigned i=0; i<data.lines.rows();++i)
    {
      lines_V_vbo.row(2*i+0) = data.lines.block<1, 3>(i, 0).cast<float>();
      lines_V_vbo.row(2*i+1) = data.lines.block<1, 3>(i, 3).cast<float>();
      lines_V_colors_vbo.row(2*i+0) = data.lines.block<1, 3>(i, 6).cast<float>();
      lines_V_colors_vbo.row(2*i+1) = data.lines.block<1, 3>(i, 6).cast<float>();
      lines_F_vbo(2*i+0) = 2*i+0;
      lines_F_vbo(2*i+1) = 2*i+1;
    }
  }

  if (dirty & ViewerData::DIRTY_OVERLAY_POINTS)
  {
    points_V_vbo.resize(data.points.rows(),3);
    points_V_colors_vbo.resize(data.points.rows(),3);
    points_F_vbo.resize(data.points.rows(),1);
    for (unsigned i=0; i<data.points.rows();++i)
    {
      points_V_vbo.row(i) = data.points.block<1, 3>(i, 0).cast<float>();
      points_V_colors_vbo.row(i) = data.points.block<1, 3>(i, 3).cast<float>();
      points_F_vbo(i) = i;
    }
  }
}

IGL_INLINE void igl::opengl::MeshGL::bind_mesh()
{
  glBindVertexArray(vao_mesh);
  glUseProgram(shader_mesh);
  bind_vertex_attrib_array(shader_mesh,"position", vbo_V, V_vbo, dirty & ViewerData::DIRTY_POSITION);
  bind_vertex_attrib_array(shader_mesh,"normal", vbo_V_normals, V_normals_vbo, dirty & ViewerData::DIRTY_NORMAL);
  bind_vertex_attrib_array(shader_mesh,"Ka", vbo_V_ambient, V_ambient_vbo, dirty & ViewerData::DIRTY_AMBIENT);
  bind_vertex_attrib_array(shader_mesh,"Kd", vbo_V_diffuse, V_diffuse_vbo, dirty & ViewerData::DIRTY_DIFFUSE);
  bind_vertex_attrib_array(shader_mesh,"Ks", vbo_V_specular, V_specular_vbo, dirty & ViewerData::DIRTY_SPECULAR);
  bind_vertex_attrib_array(shader_mesh,"texcoord", vbo_V_uv, V_uv_vbo, dirty & ViewerData::DIRTY_UV);

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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex_u, tex_v, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.data());
  }
  glUniform1i(glGetUniformLocation(shader_mesh,"tex"), 0);
  dirty &= ~ViewerData::DIRTY_MESH;
}

IGL_INLINE void igl::opengl::MeshGL::bind_overlay_lines()
{
  bool is_dirty = dirty & ViewerData::DIRTY_OVERLAY_LINES;

  glBindVertexArray(vao_overlay_lines);
  glUseProgram(shader_overlay_lines);
 bind_vertex_attrib_array(shader_overlay_lines,"position", vbo_lines_V, lines_V_vbo, is_dirty);
 bind_vertex_attrib_array(shader_overlay_lines,"color", vbo_lines_V_colors, lines_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_lines_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*lines_F_vbo.size(), lines_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~ViewerData::DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::opengl::MeshGL::bind_overlay_points()
{
  bool is_dirty = dirty & ViewerData::DIRTY_OVERLAY_POINTS;

  glBindVertexArray(vao_overlay_points);
  glUseProgram(shader_overlay_points);
 bind_vertex_attrib_array(shader_overlay_points,"position", vbo_points_V, points_V_vbo, is_dirty);
 bind_vertex_attrib_array(shader_overlay_points,"color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*points_F_vbo.size(), points_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~ViewerData::DIRTY_OVERLAY_POINTS;
}

IGL_INLINE void igl::opengl::MeshGL::draw_mesh(bool solid)
{
  glPolygonMode(GL_FRONT_AND_BACK, solid ? GL_FILL : GL_LINE);

  /* Avoid Z-buffer fighting between filled triangles & wireframe lines */
  if (solid)
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
  }
  glDrawElements(GL_TRIANGLES, 3*F_vbo.rows(), GL_UNSIGNED_INT, 0);

  glDisable(GL_POLYGON_OFFSET_FILL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

IGL_INLINE void igl::opengl::MeshGL::draw_overlay_lines()
{
  glDrawElements(GL_LINES, lines_F_vbo.rows(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::opengl::MeshGL::draw_overlay_points()
{
  glDrawElements(GL_POINTS, points_F_vbo.rows(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::opengl::MeshGL::init()
{
  if(is_initialized)
  {
    return;
  }
  is_initialized = true;
  std::string mesh_vertex_shader_string =
  "#version 150\n"
  "uniform mat4 model;"
  "uniform mat4 view;"
  "uniform mat4 proj;"
  "in vec3 position;"
  "in vec3 normal;"
  "out vec3 position_eye;"
  "out vec3 normal_eye;"
  "in vec4 Ka;"
  "in vec4 Kd;"
  "in vec4 Ks;"
  "in vec2 texcoord;"
  "out vec2 texcoordi;"
  "out vec4 Kai;"
  "out vec4 Kdi;"
  "out vec4 Ksi;"

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
  "in vec4 Ksi;"
  "in vec4 Kdi;"
  "in vec4 Kai;"
  "in vec2 texcoordi;"
  "uniform sampler2D tex;"
  "uniform float specular_exponent;"
  "uniform float lighting_factor;"
  "uniform float texture_factor;"
  "out vec4 outColor;"
  "void main()"
  "{"
  "vec3 Ia = La * vec3(Kai);"    // ambient intensity

  "vec3 light_position_eye = vec3 (view * vec4 (light_position_world, 1.0));"
  "vec3 vector_to_light_eye = light_position_eye - position_eye;"
  "vec3 direction_to_light_eye = normalize (vector_to_light_eye);"
  "float dot_prod = dot (direction_to_light_eye, normal_eye);"
  "float clamped_dot_prod = max (dot_prod, 0.0);"
  "vec3 Id = Ld * vec3(Kdi) * clamped_dot_prod;"    // Diffuse intensity

  "vec3 reflection_eye = reflect (-direction_to_light_eye, normal_eye);"
  "vec3 surface_to_viewer_eye = normalize (-position_eye);"
  "float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);"
  "dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);"
  "float specular_factor = pow (dot_prod_specular, specular_exponent);"
  "vec3 Is = Ls * vec3(Ksi) * specular_factor;"    // specular intensity
  "vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);"
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
  create_shader_program(
    mesh_vertex_shader_string,
    mesh_fragment_shader_string,
    {},
    shader_mesh);
  create_shader_program(
    overlay_vertex_shader_string,
    overlay_fragment_shader_string,
    {},
    shader_overlay_lines);
  create_shader_program(
    overlay_vertex_shader_string,
    overlay_point_fragment_shader_string,
    {},
    shader_overlay_points);
}

IGL_INLINE void igl::opengl::MeshGL::free()
{
  const auto free = [](GLuint & id)
  {
    if(id)
    {
      destroy_shader_program(id);
      id = 0;
    }
  };
  free(shader_mesh);
  free(shader_overlay_lines);
  free(shader_overlay_points);
  free_buffers();
}
