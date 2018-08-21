std::string lines_geom_shader = R"(

#version 330

layout(lines) in;
layout(triangle_strip, max_vertices=24) out;

// View parameters
uniform bool orthographic;
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
  flat vec4 a, b;       // endpoints in camera space, last coordinate = radius
  flat vec3 ex, ey, ez; // impostor frame
  flat float alpha;     // precomputed per impostor
  smooth vec4 ambient_color;
  smooth vec4 diffuse_color;
  smooth vec4 specular_color;
  smooth vec3 mapping;
} frag;

// Emit a corner vertex
float prepare()
{
  frag.a = vec4(vert[0].point_position_eye, vert[0].point_radius);
  frag.b = vec4(vert[1].point_position_eye, vert[1].point_radius);

  vec3 axis_center = vec3(mix(frag.a, frag.b, 0.5));
  vec3 ray_origin;
  vec3 ray_direction;
  if (orthographic) {
    ray_origin = vec3(axis_center.xy, 0);
    ray_direction = normalize(axis_center - ray_origin);
  } else {
    ray_origin = vec3(0.0);
    ray_direction = normalize(axis_center);
  }

  frag.ey = normalize(vec3(frag.b - frag.a));
  float sign = (dot(frag.ey, ray_direction) > 0 ? 1 : -1);
  frag.ey *= sign;
  frag.ex = normalize(cross(frag.ey, ray_direction));
  frag.ez = normalize(cross(frag.ex, frag.ey));

  // if (sign < 0) {
  //   vec4 tmp = frag.b;
  //   frag.b = frag.a;
  //   frag.a = tmp;
  // }

  float ra = frag.a.w;
  float rb = frag.b.w;
  frag.alpha = (rb - ra) / dot(vec3(frag.b - frag.a), frag.ey);

  return sign;
}

void corner(vec3 coord)
{
  frag.ambient_color  = mix(vert[0].ambient_color, vert[1].ambient_color, 0.5*coord.t+0.5);
  frag.diffuse_color  = mix(vert[0].diffuse_color, vert[1].diffuse_color, 0.5*coord.t+0.5);
  frag.specular_color = mix(vert[0].specular_color, vert[1].specular_color, 0.5*coord.t+0.5);

  vec4 p = mix(frag.a, frag.b, 0.5 * coord.t + 0.5);
  float r = p.w;
  float rmax = max(frag.a.w, frag.b.w);
  vec3 corner_position_eye = p.xyz + rmax * (coord.s * frag.ex + coord.p * frag.ez);
  gl_Position = proj * vec4(corner_position_eye, 1.0);

  frag.mapping = coord;

  gl_PrimitiveID = gl_PrimitiveIDIn;
  EmitVertex();
}

void main()
{
  float sign = prepare();
  corner(vec3(-1.0, -1.0, -1.0));
  corner(vec3(-1.0,  1.0, -1.0));
  corner(vec3( 1.0, -1.0, -1.0));
  corner(vec3( 1.0,  1.0, -1.0));
  EndPrimitive();
  // corner(vec3(-1.0, -1.0,  1.0));
  // corner(vec3(-1.0,  1.0,  1.0));
  // corner(vec3( 1.0, -1.0,  1.0));
  // corner(vec3( 1.0,  1.0,  1.0));
  // EndPrimitive();
  // corner(vec3(-1.0, -1.0, -1.0));
  // corner(vec3(-1.0,  1.0, -1.0));
  // corner(vec3(-1.0, -1.0,  1.0));
  // corner(vec3(-1.0,  1.0,  1.0));
  // EndPrimitive();
  // corner(vec3( 1.0, -1.0, -1.0));
  // corner(vec3( 1.0,  1.0, -1.0));
  // corner(vec3( 1.0, -1.0,  1.0));
  // corner(vec3( 1.0,  1.0,  1.0));
  // EndPrimitive();
  corner(vec3(-1.0, -sign, -1.0));
  corner(vec3(-1.0, -sign,  1.0));
  corner(vec3( 1.0, -sign, -1.0));
  corner(vec3( 1.0, -sign,  1.0));
  EndPrimitive();
  // corner(vec3(-1.0, sign, -1.0));
  // corner(vec3(-1.0, sign,  1.0));
  // corner(vec3( 1.0, sign, -1.0));
  // corner(vec3( 1.0, sign,  1.0));
  // EndPrimitive();
}

)";
