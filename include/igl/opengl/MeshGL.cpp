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
  glGenBuffers(1, &vbo_points_V_radius);
  glGenBuffers(1, &vbo_points_V_colors);

  dirty = MeshGL::DIRTY_ALL;
}

IGL_INLINE void igl::opengl::MeshGL::free_buffers()
{
  if (is_initialized)
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
    glDeleteBuffers(1, &vbo_points_V_radius);
    glDeleteBuffers(1, &vbo_points_V_colors);

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
  bind_vertex_attrib_array(shader_overlay_points,"radius", vbo_points_V_radius, points_V_radius_vbo, is_dirty);
  bind_vertex_attrib_array(shader_overlay_points,"color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
  if (is_dirty)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*points_F_vbo.size(), points_F_vbo.data(), GL_DYNAMIC_DRAW);

  dirty &= ~MeshGL::DIRTY_OVERLAY_POINTS;
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
  uniform vec3 light_position_eye;
  in vec3 position_eye;
  in vec3 normal_eye;
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


  std::string point_vertex_shader_string =
R"(#version 330

// View parameters
uniform float scaling_factor; // ‖ diag(view) ‖
uniform mat4 view;

// Input
layout(location = 0) in vec3  position;
layout(location = 1) in float radius;
layout(location = 2) in vec3  color;

// Output
out VertexData {
  vec3 cameraSpherePos;
  float sphereRadius;
  vec4 ambientColor;
  vec4 diffuseColor;
  vec4 specularColor;
} outData;

void main() {
  outData.cameraSpherePos = vec3(view * vec4(position.xyz, 1.0));
  outData.sphereRadius = radius * scaling_factor;
  outData.ambientColor = vec4(0.1 * color, 1.0);
  outData.diffuseColor = vec4(color, 1.0);
  vec3 grey = vec3(0.3);
  outData.specularColor = vec4(grey + 0.1 * (color.xyz- grey), 1.0);
})";

  std::string point_geom_shader_string =
R"(#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

// View parameters
uniform mat4 proj;

// Input
in VertexData {
  vec3 cameraSpherePos;
  float sphereRadius;
  vec4 ambientColor;
  vec4 diffuseColor;
  vec4 specularColor;
} vert[];

// Output
out FragData {
  flat vec3 cameraSpherePos;
  flat float sphereRadius;
  flat vec4 ambientColor;
  flat vec4 diffuseColor;
  flat vec4 specularColor;
  smooth vec2 mapping;
} frag;

// Constant
const float g_BoxCorrection = 1.5;

// Emit a corner vertex
void corner(vec2 coord) {
  frag.cameraSpherePos = vec3(vert[0].cameraSpherePos);
  frag.sphereRadius = vert[0].sphereRadius;
  frag.mapping = coord * g_BoxCorrection;
  frag.ambientColor = vert[0].ambientColor;
  frag.diffuseColor = vert[0].diffuseColor;
  frag.specularColor = vert[0].specularColor;
  vec4 cameraCornerPos = vec4(vert[0].cameraSpherePos, 1.0);
  cameraCornerPos.xy += vert[0].sphereRadius * coord * g_BoxCorrection;
  gl_Position = proj * cameraCornerPos;
  gl_PrimitiveID = gl_PrimitiveIDIn;
  EmitVertex();
}

void main() {
  corner(vec2(-1.0, -1.0)); // Bottom-left
  corner(vec2(-1.0, 1.0));  // Top-left
  corner(vec2(1.0, -1.0));  // Bottom-right
  corner(vec2(1.0, 1.0));   // Top-right
})";

  std::string point_fragment_shader_string =
R"(// -----------------------------------------------------------------------------
// Code adapted from:
// https://alfonse.bitbucket.io/oldtut/Illumination/Tutorial%2013.html
// -----------------------------------------------------------------------------

#version 330

// Material
uniform float specular_exponent;

// Light
uniform vec4 fixed_color;
uniform vec3 light_position_eye;
uniform float lighting_factor;

// View parameters
uniform mat4 proj;

// Input
in FragData {
  flat vec3 cameraSpherePos;
  flat float sphereRadius;
  flat vec4 ambientColor;
  flat vec4 diffuseColor;
  flat vec4 specularColor;
  smooth vec2 mapping;
};

// Output
layout(location = 0) out vec4 fragColor;

// -----------------------------------------------------------------------------

// Lighting equation
vec4 ComputeLighting(
  in vec3 light_position_eye,
  in vec3 position_eye, in vec3 normal_eye,
  in vec3 La, in vec3 Ld, in vec3 Ls,
  in vec4 Ka, in vec4 Kd, in vec4 Ks)
{
  vec3 Ia = La * vec3(Ka);    // ambient intensity

  vec3 vector_to_light_eye = light_position_eye - position_eye;
  vec3 direction_to_light_eye = normalize(vector_to_light_eye);
  float dot_prod = dot (direction_to_light_eye, normal_eye);
  float clamped_dot_prod = max (dot_prod, 0.0);
  vec3 Id = vec3(Kd) * clamped_dot_prod;    // Diffuse intensity

  vec3 reflection_eye = reflect (-direction_to_light_eye, normal_eye);
  vec3 surface_to_viewer_eye = normalize (-position_eye);
  float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);
  dot_prod_specular = float(abs(dot_prod)==dot_prod) * max (dot_prod_specular, 0.0);
  float specular_factor = pow (dot_prod_specular, specular_exponent);
  vec3 Is = Ls * vec3(Ks) * specular_factor;    // specular intensity
  vec4 outColor = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Ka), (Ka.a+Ks.a+Kd.a)/3);
  if (fixed_color != vec4(0.0)) outColor = fixed_color;
  return outColor;
}

// Compute sphere position and normal in camera space
void Impostor(out vec3 cameraPos, out vec3 cameraNormal) {
  vec3 cameraPlanePos = vec3(mapping * sphereRadius, 0.0) + cameraSpherePos;
  vec3 rayDirection = normalize(cameraPlanePos);

  float B = 2.0 * dot(rayDirection, -cameraSpherePos);
  float C = dot(cameraSpherePos, cameraSpherePos) - (sphereRadius * sphereRadius);

  float det = (B * B) - (4 * C);
  if (det < 0.0) {
    discard;
  }

  float sqrtDet = sqrt(det);
  float posT = (-B + sqrtDet)/2;
  float negT = (-B - sqrtDet)/2;

  float intersectT = min(posT, negT);
  cameraPos = rayDirection * intersectT;
  cameraNormal = normalize(cameraPos - cameraSpherePos);
}

void main() {
  vec3 cameraPos;
  vec3 cameraNormal;

  Impostor(cameraPos, cameraNormal);

  // Set the depth based on the new cameraPos.
  vec4 clipPos = proj * vec4(cameraPos, 1.0);
  float ndcDepth = clipPos.z / clipPos.w;
  gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

  // Compute lighting
  fragColor = ComputeLighting(
    light_position_eye, cameraPos, cameraNormal,
    vec3(1.0), vec3(1.0), vec3(1.0),
    ambientColor, diffuseColor, specularColor);
})";

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
    point_geom_shader_string,
    point_vertex_shader_string,
    point_fragment_shader_string,
    {},
    shader_overlay_points);
}

IGL_INLINE void igl::opengl::MeshGL::free()
{
  const auto free_shader = [](GLuint & id)
  {
    if(id)
    {
      destroy_shader_program(id);
      id = 0;
    }
  };

  if (is_initialized)
  {
    free_shader(shader_mesh);
    free_shader(shader_overlay_lines);
    free_shader(shader_overlay_points);
    free_buffers();
  }
}
