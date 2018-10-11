std::string mesh_vertex_shader = R"(

#version 330

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

out vec4 Kai;
out vec4 Kdi;
out vec4 Ksi;
out vec2 texcoordi;

void main()
{
  position_eye = vec3(view * vec4(position, 1.0));
  normal_eye = vec3(normal_matrix * vec4(normal, 0.0));
  normal_eye = normalize(normal_eye);
  gl_Position = proj * vec4(position_eye, 1.0);
  Kai = Ka;
  Kdi = Kd;
  Ksi = Ks;
  texcoordi = texcoord;
}

)";
