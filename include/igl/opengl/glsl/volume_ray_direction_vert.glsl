// Shader transforming the vertices from model coordinates to clip space.
// We use this to render two framebuffers containing the intersection of a
// ray entering and exiting a bounding box respectively

#version 150
in vec3 in_position;
out vec3 color;

uniform mat4 model_matrix;
uniform mat4 view_matrix;
uniform mat4 projection_matrix;

void main() {
  gl_Position = projection_matrix * view_matrix * model_matrix * vec4(in_position, 1.0);
  color = in_position.xyz;
}
