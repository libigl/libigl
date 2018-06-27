// Vertex shader that is used to trigger the volume rendering by rendering a static
// screen-space filling quad.

#version 150

// Create two triangles that are filling the entire screen [-1, 1]
vec2 positions[6] = vec2[](
  vec2(-1.0, -1.0),
  vec2( 1.0, -1.0),
  vec2( 1.0,  1.0),

  vec2(-1.0, -1.0),
  vec2( 1.0,  1.0),
  vec2(-1.0,  1.0)
);

out vec2 uv;

void main() {
  // Clipspace in [-1, 1]
  gl_Position = vec4(positions[gl_VertexID], 0.0, 1.0);

  // UV coordinates in [0, 1]
  uv = (positions[gl_VertexID] + 1.0) / 2.0;
}
