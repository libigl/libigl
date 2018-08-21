std::string mesh_fragment_shader = R"(

#version 330

// Camera
uniform mat4 view;
uniform mat4 proj;

// Lights
uniform vec4 fixed_color;
uniform vec3 light_vector_eye;
uniform bool orthographic;

const vec3 Ls = vec3 (1, 1, 1);
const vec3 Ld = vec3 (1, 1, 1);
const vec3 La = vec3 (1, 1, 1);

// Surface
in vec3 position_eye;
in vec3 normal_eye;

in vec4 Ksi;
in vec4 Kdi;
in vec4 Kai;
in vec2 texcoordi;

// Misc
uniform sampler2D tex;
uniform float specular_exponent;
uniform float lighting_factor;
uniform float texture_factor;
out vec4 outColor;

void main()
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

  vec4 color = vec4(lighting_factor * (Is + Id) + Ia + (1.0-lighting_factor) * vec3(Kdi),(Kai.a+Ksi.a+Kdi.a)/3);
  outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;
  if (fixed_color != vec4(0.0)) outColor = fixed_color;
}

)";
