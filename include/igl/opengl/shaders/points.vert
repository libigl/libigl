std::string points_vertex_shader = R"(

#version 330

// View parameters
uniform float scaling_factor;
uniform mat4 view;

// Input
layout(location = 0) in vec3  position;
layout(location = 1) in float radius;
layout(location = 2) in vec3  color;

// Output
out VertexData {
  vec3 sphere_position_eye;
  float sphere_radius;
  vec4 ambient_color;
  vec4 diffuse_color;
  vec4 specular_color;
} vertex;

void main()
{
  vertex.sphere_position_eye = vec3(view * vec4(position.xyz, 1.0));
  vertex.sphere_radius = radius * scaling_factor;
  vertex.ambient_color = vec4(0.1 * color, 1.0);
  vertex.diffuse_color = vec4(color, 1.0);
  vec3 grey = vec3(0.3);
  vertex.specular_color = vec4(grey + 0.1 * (color.xyz - grey), 1.0);
}

)";
