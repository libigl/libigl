std::string points_fragment_shader = R"(

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
  flat vec3 sphere_position_eye;
  flat float sphere_radius;
  flat vec4 ambient_color;
  flat vec4 diffuse_color;
  flat vec4 specular_color;
  smooth vec2 mapping;
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

// Compute sphere position and normal in camera space
void impostor(out vec3 hit_position_eye, out vec3 hit_normal_eye)
{
  vec3 hit_plane_position_eye = vec3(mapping * sphere_radius, 0.0) + sphere_position_eye;

  vec3 ray_origin;
  vec3 ray_direction;
  if (orthographic)
  {
    ray_origin = vec3(hit_plane_position_eye.xy, 0);
    ray_direction = normalize(hit_plane_position_eye - ray_origin);
  }
  else
  {
    ray_origin = vec3(0.0);
    ray_direction = normalize(hit_plane_position_eye);
  }
  vec3 camera_to_sphere = ray_origin - sphere_position_eye;

  float B = 2.0 * dot(ray_direction, camera_to_sphere);
  float C = dot(camera_to_sphere, camera_to_sphere) - (sphere_radius * sphere_radius);

  float det = (B * B) - (4 * C);
  if (det < 0.0)
  {
    discard;
  }

  float sqrt_det = sqrt(det);
  float t1 = (-B + sqrt_det)/2;
  float t2 = (-B - sqrt_det)/2;
  float t = min(t1, t2);

  hit_position_eye = ray_origin + ray_direction * t;
  hit_normal_eye = normalize(hit_position_eye - sphere_position_eye);
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

  // Compute lighting
  frag_color = compute_lighting(
    hit_position_eye, hit_normal_eye,
    vec3(1.0), vec3(1.0), vec3(1.0),
    ambient_color, diffuse_color, specular_color);
}

)";
