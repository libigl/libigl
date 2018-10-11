std::string points_geom_shader = R"(

#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

// View parameters
uniform mat4 proj;

// Input
in VertexData {
  vec3 sphere_position_eye;
  float sphere_radius;
  vec4 ambient_color;
  vec4 diffuse_color;
  vec4 specular_color;
} vert[];

// Output
out FragData {
  flat vec3 sphere_position_eye;
  flat float sphere_radius;
  flat vec4 ambient_color;
  flat vec4 diffuse_color;
  flat vec4 specular_color;
  smooth vec2 mapping;
} frag;

// Constant
const float BOX_CORRECTION = 1.5;

// Emit a corner vertex
void corner(vec2 coord)
{
  frag.sphere_position_eye = vec3(vert[0].sphere_position_eye);
  frag.sphere_radius = vert[0].sphere_radius;
  frag.mapping = coord * BOX_CORRECTION;
  frag.ambient_color = vert[0].ambient_color;
  frag.diffuse_color = vert[0].diffuse_color;
  frag.specular_color = vert[0].specular_color;
  vec4 corner_position_eye = vec4(vert[0].sphere_position_eye, 1.0);
  corner_position_eye.xy += vert[0].sphere_radius * coord * BOX_CORRECTION;
  gl_Position = proj * corner_position_eye;
  gl_PrimitiveID = gl_PrimitiveIDIn;
  EmitVertex();
}

void main()
{
  corner(vec2(-1.0, -1.0)); // Bottom-left
  corner(vec2(-1.0, 1.0));  // Top-left
  corner(vec2(1.0, -1.0));  // Bottom-right
  corner(vec2(1.0, 1.0));   // Top-right
}

)";
