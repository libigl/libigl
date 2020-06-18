// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "MeshGL.h"
#include "bind_vertex_attrib_array.h"
#include "create_shader_program.h"
#include "destroy_shader_program.h"
#include "shaders/text.vert"
#include "shaders/text.geom"
#include "shaders/text.frag"
#include <igl/png/texture_from_png.h>
#include <iostream>

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

  // Vert ID Labels
  glGenVertexArrays(1, &vao_vid_labels);
  glBindVertexArray(vao_vid_labels);
  glGenBuffers(1, &vbo_vid_labels_pos);
  glGenBuffers(1, &vbo_vid_labels_characters);
  glGenBuffers(1, &vbo_vid_labels_offset);
  glGenBuffers(1, &vbo_vid_labels_indices);
  // Face ID Labels
  glGenVertexArrays(1, &vao_fid_labels);
  glBindVertexArray(vao_fid_labels);
  glGenBuffers(1, &vbo_fid_labels_pos);
  glGenBuffers(1, &vbo_fid_labels_characters);
  glGenBuffers(1, &vbo_fid_labels_offset);
  glGenBuffers(1, &vbo_fid_labels_indices);
  glGenVertexArrays(1, &vao_extra_labels);
  // Extra Labels
  glBindVertexArray(vao_extra_labels);
  glGenBuffers(1, &vbo_extra_labels_pos);
  glGenBuffers(1, &vbo_extra_labels_characters);
  glGenBuffers(1, &vbo_extra_labels_offset);
  glGenBuffers(1, &vbo_extra_labels_indices);

  dirty = MeshGL::DIRTY_ALL;
}

IGL_INLINE void igl::opengl::MeshGL::free_buffers()
{
  if (is_initialized)
  {
    glDeleteVertexArrays(1, &vao_mesh);
    glDeleteVertexArrays(1, &vao_overlay_lines);
    glDeleteVertexArrays(1, &vao_overlay_points);
    glDeleteVertexArrays(1, &vao_vid_labels);
    glDeleteVertexArrays(1, &vao_fid_labels);
    glDeleteVertexArrays(1, &vao_extra_labels);

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
    
    // Vert ID Labels
    glDeleteBuffers(1, &vbo_vid_labels_pos);
    glDeleteBuffers(1, &vbo_vid_labels_characters);
    glDeleteBuffers(1, &vbo_vid_labels_offset);
    glDeleteBuffers(1, &vbo_vid_labels_indices);
    // Face ID Labels
    glDeleteBuffers(1, &vbo_fid_labels_pos);
    glDeleteBuffers(1, &vbo_fid_labels_characters);
    glDeleteBuffers(1, &vbo_fid_labels_offset);
    glDeleteBuffers(1, &vbo_fid_labels_indices);
    // Extra Labels 
    glDeleteBuffers(1, &vbo_extra_labels_pos);
    glDeleteBuffers(1, &vbo_extra_labels_characters);
    glDeleteBuffers(1, &vbo_extra_labels_offset);
    glDeleteBuffers(1, &vbo_extra_labels_indices);

    glDeleteTextures(1, &vbo_tex);
  }
}

IGL_INLINE void igl::opengl::MeshGL::bind_mesh()
{
  glBindVertexArray(vao_mesh);
  glUseProgram(shader_mesh);
  bind_vertex_attrib_array(shader_mesh,"position", vbo_V, V_vbo, dirty & MeshGL::DIRTY_POSITION);
  bind_vertex_attrib_array(shader_mesh,"normal", vbo_V_normals, V_normals_vbo, dirty & MeshGL::DIRTY_NORMAL);
  bind_vertex_attrib_array(shader_mesh,"Ka", vbo_V_ambient, V_ambient_vbo, dirty & MeshGL::DIRTY_AMBIENT);
  bind_vertex_attrib_array(shader_mesh,"Kd", vbo_V_diffuse, V_diffuse_vbo, dirty & MeshGL::DIRTY_DIFFUSE);
  bind_vertex_attrib_array(shader_mesh,"Ks", vbo_V_specular, V_specular_vbo, dirty & MeshGL::DIRTY_SPECULAR);
  bind_vertex_attrib_array(shader_mesh,"texcoord", vbo_V_uv, V_uv_vbo, dirty & MeshGL::DIRTY_UV);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_F);
  if (dirty & MeshGL::DIRTY_FACE)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*F_vbo.size(), F_vbo.data(), GL_DYNAMIC_DRAW);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, vbo_tex);
  if (dirty & MeshGL::DIRTY_TEXTURE)
  {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex_u, tex_v, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.data());
  }
  glUniform1i(glGetUniformLocation(shader_mesh,"tex"), 0);
  dirty &= ~MeshGL::DIRTY_MESH;
}

IGL_INLINE void igl::opengl::MeshGL::bind_overlay_lines()
{
  bool is_dirty = dirty & MeshGL::DIRTY_OVERLAY_LINES;

  glBindVertexArray(vao_overlay_lines);
  glUseProgram(shader_overlay_lines);
 bind_vertex_attrib_array(shader_overlay_lines,"position", vbo_lines_V, lines_V_vbo, is_dirty);
 bind_vertex_attrib_array(shader_overlay_lines,"color", vbo_lines_V_colors, lines_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_lines_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*lines_F_vbo.size(), lines_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~MeshGL::DIRTY_OVERLAY_LINES;
}

IGL_INLINE void igl::opengl::MeshGL::bind_overlay_points()
{
  bool is_dirty = dirty & MeshGL::DIRTY_OVERLAY_POINTS;

  glBindVertexArray(vao_overlay_points);
  glUseProgram(shader_overlay_points);
 bind_vertex_attrib_array(shader_overlay_points,"position", vbo_points_V, points_V_vbo, is_dirty);
 bind_vertex_attrib_array(shader_overlay_points,"color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*points_F_vbo.size(), points_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~MeshGL::DIRTY_OVERLAY_POINTS;
}

IGL_INLINE void igl::opengl::MeshGL::bind_font_atlas()
{
  const std::string font_atlas = "/home/michelle/Documents/LIBIGL/opengl_text_rendering/libigl/include/igl/opengl/shaders/verasansmono.png";
  GLuint texture_handle;
  igl::png::texture_from_png(font_atlas, texture_handle);
  glBindTexture(GL_TEXTURE_2D, texture_handle);
}

IGL_INLINE void igl::opengl::MeshGL::bind_vid_labels()
{
  bool is_dirty = dirty & MeshGL::DIRTY_VID_LABELS;
  glBindVertexArray(vao_vid_labels);
  glUseProgram(shader_text);
  bind_vertex_attrib_array(shader_text, "position" , vbo_vid_labels_pos       , vid_label_pos_vbo   , is_dirty);
  bind_vertex_attrib_array(shader_text, "character", vbo_vid_labels_characters, vid_label_char_vbo  , is_dirty);
  bind_vertex_attrib_array(shader_text, "offset"   , vbo_vid_labels_offset    , vid_label_offset_vbo, is_dirty);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_vid_labels_indices);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*vid_label_indices_vbo.size(), vid_label_indices_vbo.data(), GL_DYNAMIC_DRAW);
  dirty &= ~MeshGL::DIRTY_VID_LABELS;
}

IGL_INLINE void igl::opengl::MeshGL::bind_fid_labels()
{
  bool is_dirty = dirty & MeshGL::DIRTY_FID_LABELS;
  glBindVertexArray(vao_fid_labels);
  glUseProgram(shader_text);
  bind_vertex_attrib_array(shader_text, "position" , vbo_fid_labels_pos       , fid_label_pos_vbo   , is_dirty);
  bind_vertex_attrib_array(shader_text, "character", vbo_fid_labels_characters, fid_label_char_vbo  , is_dirty);
  bind_vertex_attrib_array(shader_text, "offset"   , vbo_fid_labels_offset    , fid_label_offset_vbo, is_dirty);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_fid_labels_indices);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*fid_label_indices_vbo.size(), fid_label_indices_vbo.data(), GL_DYNAMIC_DRAW);
  dirty &= ~MeshGL::DIRTY_FID_LABELS;
}

IGL_INLINE void igl::opengl::MeshGL::bind_extra_labels()
{
  bool is_dirty = dirty & MeshGL::DIRTY_EXTRA_LABELS;
  glBindVertexArray(vao_extra_labels);
  glUseProgram(shader_text);
  bind_vertex_attrib_array(shader_text, "position" , vbo_extra_labels_pos       , extra_label_pos_vbo   , is_dirty);
  bind_vertex_attrib_array(shader_text, "character", vbo_extra_labels_characters, extra_label_char_vbo  , is_dirty);
  bind_vertex_attrib_array(shader_text, "offset"   , vbo_extra_labels_offset    , extra_label_offset_vbo, is_dirty);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_extra_labels_indices);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*extra_label_indices_vbo.size(), extra_label_indices_vbo.data(), GL_DYNAMIC_DRAW);
  dirty &= ~MeshGL::DIRTY_EXTRA_LABELS;
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

IGL_INLINE void igl::opengl::MeshGL::draw_vid_labels()
{
  glDrawElements(GL_POINTS, vid_label_indices_vbo.rows(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::opengl::MeshGL::draw_fid_labels()
{
  glDrawElements(GL_POINTS, fid_label_indices_vbo.rows(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::opengl::MeshGL::draw_extra_labels()
{
  glDrawElements(GL_POINTS, extra_label_indices_vbo.rows(), GL_UNSIGNED_INT, 0);
}

IGL_INLINE void igl::opengl::MeshGL::init()
{
  if(is_initialized)
  {
    return;
  }
  is_initialized = true;
  std::string mesh_vertex_shader_string =
R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform mat4 normal_matrix;
  in vec3 position;
  in vec3 normal;
  out vec3 position_eye;
  out vec3 normal_eye;
  in vec4 Ka;
  in vec4 Kd;
  in vec4 Ks;
  in vec2 texcoord;
  out vec2 texcoordi;
  out vec4 Kai;
  out vec4 Kdi;
  out vec4 Ksi;

  void main()
  {
    position_eye = vec3 (view * vec4 (position, 1.0));
    normal_eye = vec3 (normal_matrix * vec4 (normal, 0.0));
    normal_eye = normalize(normal_eye);
    gl_Position = proj * vec4 (position_eye, 1.0); //proj * view * vec4(position, 1.0);"
    Kai = Ka;
    Kdi = Kd;
    Ksi = Ks;
    texcoordi = texcoord;
  }
)";

  std::string mesh_fragment_shader_string =
R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  uniform vec4 fixed_color;
  in vec3 position_eye;
  in vec3 normal_eye;
  uniform vec3 light_position_eye;
  vec3 Ls = vec3 (1, 1, 1);
  vec3 Ld = vec3 (1, 1, 1);
  vec3 La = vec3 (1, 1, 1);
  in vec4 Ksi;
  in vec4 Kdi;
  in vec4 Kai;
  in vec2 texcoordi;
  uniform sampler2D tex;
  uniform float specular_exponent;
  uniform float lighting_factor;
  uniform float texture_factor;
  out vec4 outColor;
  void main()
  {
    vec3 Ia = La * vec3(Kai);    // ambient intensity

    vec3 vector_to_light_eye = light_position_eye - position_eye;
    vec3 direction_to_light_eye = normalize (vector_to_light_eye);
    float dot_prod = dot (direction_to_light_eye, normalize(normal_eye));
    float clamped_dot_prod = max (dot_prod, 0.0);
    vec3 Id = Ld * vec3(Kdi) * clamped_dot_prod;    // Diffuse intensity

    vec3 reflection_eye = reflect (-direction_to_light_eye, normalize(normal_eye));
    vec3 surface_to_viewer_eye = normalize (-position_eye);
    float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
    dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);
    float specular_factor = pow (dot_prod_specular, specular_exponent);
    vec3 Is = Ls * vec3(Ksi) * specular_factor;    // specular intensity
    vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);
    outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;
    if (fixed_color != vec4(0.0)) outColor = fixed_color;
  }
)";

  std::string overlay_vertex_shader_string =
R"(#version 150
  uniform mat4 view;
  uniform mat4 proj;
  in vec3 position;
  in vec3 color;
  out vec3 color_frag;

  void main()
  {
    gl_Position = proj * view * vec4 (position, 1.0);
    color_frag = color;
  }
)";

  std::string overlay_fragment_shader_string =
R"(#version 150
  in vec3 color_frag;
  out vec4 outColor;
  void main()
  {
    outColor = vec4(color_frag, 1.0);
  }
)";

  std::string overlay_point_fragment_shader_string =
R"(#version 150
  in vec3 color_frag;
  out vec4 outColor;
  void main()
  {
    if (length(gl_PointCoord - vec2(0.5)) > 0.5)
      discard;
    outColor = vec4(color_frag, 1.0);
  }
)";

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
  create_shader_program(
    text_geom_shader,
    text_vert_shader,
    text_frag_shader,
    {},
    shader_text);
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

  if (is_initialized)
  {
    free(shader_mesh);
    free(shader_overlay_lines);
    free(shader_overlay_points);
    free_buffers();
  }
}
