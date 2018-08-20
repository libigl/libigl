std::string lines_geom_shader = R"(

#version 330

layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

// View parameters
uniform mat4 proj;

// Input
in VertexData {
  vec3 point_position_eye;
  float point_radius;
  vec4 ambient_color;
  vec4 diffuse_color;
  vec4 specular_color;
} vert[];

// Output
out FragData {
  flat vec4 a, b;         // endpoints in camera space, last coordinate = radius
  flat vec3 ex, ey, ez;   // impostor frame
  flat float alpha, beta; // precomputed per impostor
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
  frag.a = vec4(vert[0].point_position_eye, vert[0].point_radius);
  frag.b = vec4(vert[1].point_position_eye, vert[1].point_radius);
  frag.mapping = coord * BOX_CORRECTION;

  frag.ey = normalize(vert[1].point_position_eye - vert[0].point_position_eye);
  frag.ez = normalize(cross(frag.ey, vec3(-frag.ey.y, frag.ey.x, 0)));
  frag.ex = normalize(cross(frag.ey, frag.ez));

  float ra = frag.a.w;
  float rb = frag.b.w;
  float ay = dot(vec3(frag.a), frag.ey);
  float by = dot(vec3(frag.b), frag.ey);
  frag.alpha = (rb - ra) / (by - ay);
  frag.beta = (ra * by - ay * rb) / (by - ay);

  frag.ambient_color  = vert[0].ambient_color; //mix(vert[0].ambient_color, vert[1].ambient_color, 0.5*coord.t+0.5);
  frag.diffuse_color  = vert[0].diffuse_color; //mix(vert[0].diffuse_color, vert[1].diffuse_color, 0.5*coord.t+0.5);
  frag.specular_color = vert[0].specular_color; //mix(vert[0].specular_color, vert[1].specular_color, 0.5*coord.t+0.5);

  vec3 corner_position_eye = mix(vert[0].point_position_eye, vert[1].point_position_eye, 0.5*coord.t*BOX_CORRECTION+0.5);
  corner_position_eye += coord.s * BOX_CORRECTION * max(vert[0].point_radius, vert[1].point_radius) * frag.ex;
  gl_Position = proj * vec4(corner_position_eye, 1.0);

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
