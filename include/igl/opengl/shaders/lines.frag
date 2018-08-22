std::string lines_fragment_shader = R"(

// -----------------------------------------------------------------------------
// Code adapted from:
// https://alfonse.bitbucket.io/oldtut/Illumination/Tutorial%2013.html
// -----------------------------------------------------------------------------

#version 330

// Material
uniform float specular_exponent;

// Light
uniform vec4 fixed_color;
uniform vec3 light_vector_eye;
uniform float lighting_factor;

// View parameters
uniform bool orthographic;
uniform mat4 proj;

// Input
in FragData {
  flat vec4 a, b;         // endpoints in camera space, last coordinate = radius
  flat vec3 ex, ey, ez;   // impostor frame
  flat float alpha;       // precomputed per impostor
  smooth vec4 ambient_color;
  smooth vec4 diffuse_color;
  smooth vec4 specular_color;
  smooth vec3 mapping;
};

// Output
layout(location = 0) out vec4 frag_color;

// -----------------------------------------------------------------------------

// Lighting equation
vec4 compute_lighting(
  in vec3 position_eye, in vec3 normal_eye,
  in vec3 La, in vec3 Ld, in vec3 Ls,
  in vec4 Kai, in vec4 Kdi, in vec4 Ksi)
{
  vec3 Ia = La * vec3(Kai); // ambient intensity

  vec3 direction_to_light_eye = (orthographic ? -light_vector_eye : normalize(light_vector_eye - position_eye));
  float dot_prod = dot(direction_to_light_eye, normalize(normal_eye));
  float clamped_dot_prod = clamp(dot_prod, 0, 1);
  vec3 Id = Ld * vec3(Kdi) * clamped_dot_prod; // diffuse intensity

  vec3 reflection_eye = reflect(-direction_to_light_eye, normalize(normal_eye));
  vec3 surface_to_viewer_eye = (orthographic ? vec3(0,0,1) : normalize(-position_eye));
  float dot_prod_specular = dot(reflection_eye, surface_to_viewer_eye);
  dot_prod_specular = float(abs(dot_prod)==dot_prod) * max(dot_prod_specular, 0.0);
  float specular_factor = pow(dot_prod_specular, specular_exponent);
  vec3 Is = Ls * vec3(Ksi) * specular_factor; // specular intensity

  vec4 outColor = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);
  if (fixed_color != vec4(0.0)) outColor = fixed_color;
  return outColor;
}

// Compute cylinder position and normal in camera space
void impostor(out vec3 hit_position_eye, out vec3 hit_normal_eye)
{
  vec4 p = mix(a, b, 0.5*mapping.t+0.5);
  float radius = p.w;
  vec3 hit_plane_position_eye = p.xyz + max(a.w, b.w) * (mapping.s * ex + mapping.p * ez);
  vec3 midpoint = vec3(mix(a, b, 0.5));
  float half_len = abs(dot(midpoint - a.xyz, ey));

  vec3 ray_origin;
  vec3 ray_direction;
  if (orthographic) {
    ray_origin = vec3(hit_plane_position_eye.xy, 0);
    ray_direction = normalize(hit_plane_position_eye - ray_origin);
  } else {
    ray_origin = vec3(0.0);
    ray_direction = normalize(hit_plane_position_eye);
  }

  vec3 ctr = vec3(dot(a.xyz, ex), dot(a.xyz, ey), dot(a.xyz, ez));
  vec3 org = vec3(dot(ray_origin, ex), dot(ray_origin, ey), dot(ray_origin, ez));
  vec3 dir = vec3(dot(ray_direction, ex), dot(ray_direction, ey), dot(ray_direction, ez));

  float A = dot(dir.xz, dir.xz) - dot(dir.y*alpha, dir.y*alpha);
  float B = 2.0 * (dot(org.xz - ctr.xz, dir.xz) - dir.y*alpha * (a.w + (org.y - ctr.y)*alpha));
  float C = dot(org.xz - ctr.xz, org.xz - ctr.xz) - dot(a.w + (org.y - ctr.y)*alpha, a.w + (org.y - ctr.y)*alpha);
  // float A = dot(dir.xz, dir.xz);
  // float B = 2.0 * dot(org.xz - ctr.xz, dir.xz);
  // float C = dot(org.xz - ctr.xz, org.xz - ctr.xz) - radius*radius;

  float det = (B * B) - (4 * A * C);
  if (det < 0.0) {
    discard;
  }

  float sqrt_det = sqrt(det);
  float t1 = (-B + sqrt_det)/(2.0 * A);
  float t2 = (-B - sqrt_det)/(2.0 * A);

  vec3 h1 = ray_origin + ray_direction * t1;
  vec3 h2 = ray_origin + ray_direction * t2;

  float t;
  if (abs(dot(h1 - midpoint, ey)) > half_len) {
    t = t2;
  } else if (abs(dot(h2 - midpoint, ey)) > half_len) {
    t = t1;
  } else {
    t = (min(t1, t2) < 0 ? max(t1, t2) : min(t1, t2));
  }

  hit_position_eye = ray_origin + ray_direction * t;

  vec3 axis_to_hit = hit_position_eye - midpoint;
  axis_to_hit = vec3(dot(axis_to_hit, ex), dot(axis_to_hit, ey), dot(axis_to_hit, ez));
  hit_normal_eye = normalize(axis_to_hit.x * ex + axis_to_hit.z * ez);
  vec3 tangent = (b.xyz + b.w  * hit_normal_eye) - (a.xyz + a.w  * hit_normal_eye);
  vec3 normal = normalize(cross(normalize(cross(tangent, hit_normal_eye)), tangent));
  hit_normal_eye = dot(normal, hit_normal_eye) > 0 ? normal : -normal;

  float t3 = dot(hit_plane_position_eye - ray_origin, ray_direction);
  vec2 axis_to_hit_plane = vec2(dot(hit_plane_position_eye - midpoint, ex), dot(hit_plane_position_eye - midpoint, ez));
  float cap_sq_radius = dot(axis_to_hit_plane, axis_to_hit_plane);
  if (abs(axis_to_hit.y) > half_len || (abs(mapping.t) == 1.0 && t3 < t && cap_sq_radius < radius*radius)) {
    if (cap_sq_radius > radius * radius) {
      discard;
    }
    hit_position_eye = hit_plane_position_eye;
    hit_normal_eye = -ey;
  }
}

void main()
{
  vec3 hit_position_eye;
  vec3 hit_normal_eye;

  impostor(hit_position_eye, hit_normal_eye);

  // Set the depth based on the new hit_position_eye.
  vec4 clipPos = proj * vec4(hit_position_eye, 1.0);
  float ndcDepth = clipPos.z / clipPos.w;
  gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

  // frag_color = vec4(abs(hit_normal_eye), 1);
  // return;

  // Compute lighting
  frag_color = compute_lighting(
    hit_position_eye, hit_normal_eye,
    vec3(1.0), vec3(1.0), vec3(1.0),
    ambient_color, diffuse_color, specular_color);
}

)";
