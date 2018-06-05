#include "VolumeGL.h"
#include "ViewerCore.h"

#include <GLFW/glfw3.h>
#include <iostream>

using namespace std;


GLuint igl::opengl::VolumeGL::BoundingBox::program = 0;
igl::opengl::VolumeGL::BoundingBox::UniformLocation igl::opengl::VolumeGL::BoundingBox::uniform_location;
GLuint igl::opengl::VolumeGL::VolumeRendering::program = 0;
igl::opengl::VolumeGL::VolumeRendering::Uniform_Location igl::opengl::VolumeGL::VolumeRendering::uniform_location;

void igl::opengl::VolumeGL::free() {
  if (!_is_initialized) {
    return;
  }
  // Bail out if there's no GL context
  if (nullptr == glfwGetCurrentContext()) {
    return;
  }

  glDeleteVertexArrays(1, &bounding_box.vao);
  glDeleteBuffers(1, &bounding_box.vbo);
  glDeleteBuffers(1, &bounding_box.ibo);
  glDeleteTextures(1, &bounding_box.entry_texture);
  glDeleteTextures(1, &bounding_box.exit_texture);
  glDeleteFramebuffers(1, &bounding_box.entry_framebuffer);
  glDeleteFramebuffers(1, &bounding_box.exit_framebuffer);

  glDeleteTextures(1, &volume_rendering.volume_texture);
  glDeleteTextures(1, &volume_rendering.transfer_function_texture);
  _is_initialized = false;
}

// TODO: Make resizing actually work
// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::resize_framebuffer_textures(const Eigen::RowVector4f& viewport) {
  cout << "resize framebuffer textures..." << endl;
  cout << "viewport is " << viewport << endl;
  glBindTexture(GL_TEXTURE_2D, bounding_box.entry_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, viewport[2],
               viewport[3], 0, GL_RGBA, GL_FLOAT, nullptr);

  glBindTexture(GL_TEXTURE_2D, bounding_box.entry_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, viewport[2],            viewport[3], 0, GL_RGBA, GL_FLOAT, nullptr);

  glBindTexture(GL_TEXTURE_2D, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::init(igl::opengl::ViewerCore& core) {
  cout << "init called!" << endl;
  // This should be enabled by default
  glEnable(GL_TEXTURE_1D);
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_TEXTURE_3D);

  //
  //   Bounding box information
  //
  glGenVertexArrays(1, &bounding_box.vao);
  glBindVertexArray(bounding_box.vao);

  // Creating the vertex buffer object
  glGenBuffers(1, &bounding_box.vbo);
  glBindBuffer(GL_ARRAY_BUFFER, bounding_box.vbo);

  // Unit cube centered around 0.5 \in [0,1]
  const GLfloat vertexData[] = {
      0.f, 0.f, 0.f,
      0.f, 0.f, 1.f,
      0.f, 1.f, 0.f,
      0.f, 1.f, 1.f,
      1.f, 0.f, 0.f,
      1.f, 0.f, 1.f,
      1.f, 1.f, 0.f,
      1.f, 1.f, 1.f
  };
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, false, 3 * sizeof(GLfloat), nullptr);
  glEnableVertexAttribArray(0);

  // Create the index buffer object for the bounding box
  glGenBuffers(1, &bounding_box.ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bounding_box.ibo);

  // Specify the 12 faces of the unit cube
  const GLubyte iboData[] = {
      0, 6, 4,
      0, 2, 6,
      0, 3, 2,
      0, 1, 3,
      2, 7, 6,
      2, 3, 7,
      4, 6, 7,
      4, 7, 5,
      0, 4, 5,
      0, 5, 1,
      1, 5, 7,
      1, 7, 3
  };
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(iboData), iboData, GL_STATIC_DRAW);
  glBindVertexArray(0);

  // Entry point texture and frame buffer
  glGenTextures(1, &bounding_box.entry_texture);
  glBindTexture(GL_TEXTURE_2D, bounding_box.entry_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, core.viewport[2],
               core.viewport[3], 0, GL_RGBA, GL_FLOAT, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glGenFramebuffers(1, &bounding_box.entry_framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, bounding_box.entry_framebuffer);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         bounding_box.entry_texture, 0);

  // Exit point texture and frame buffer
  glGenTextures(1, &bounding_box.exit_texture);
  glBindTexture(GL_TEXTURE_2D, bounding_box.exit_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, core.viewport[2],
               core.viewport[3], 0, GL_RGBA, GL_FLOAT, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glGenFramebuffers(1, &bounding_box.exit_framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, bounding_box.exit_framebuffer);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         bounding_box.exit_texture, 0);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);

  // Volume texture
  glGenTextures(1, &volume_rendering.volume_texture);
  glBindTexture(GL_TEXTURE_3D, volume_rendering.volume_texture);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);


  // Transfer function texture
  glGenTextures(1, &volume_rendering.transfer_function_texture);
  glBindTexture(GL_TEXTURE_1D, volume_rendering.transfer_function_texture);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

  _is_initialized = true;
  return true;
}


// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::set_data(const Eigen::RowVector3i& dimensions, const Eigen::VectorXd& data) {
  volume_rendering_parameters.volume_dimensions = {
    static_cast<GLuint>(dimensions[0]),
    static_cast<GLuint>(dimensions[1]),
    static_cast<GLuint>(dimensions[2]) };
  volume_rendering_parameters.volume_dimensions_rcp = {
      1.f / volume_rendering_parameters.volume_dimensions[0],
      1.f / volume_rendering_parameters.volume_dimensions[1],
      1.f / volume_rendering_parameters.volume_dimensions[2]
  };
  GLuint maxDim = *std::max_element(
      volume_rendering_parameters.volume_dimensions.begin(),
      volume_rendering_parameters.volume_dimensions.end());
  volume_rendering_parameters.normalized_volume_dimensions = {
      volume_rendering_parameters.volume_dimensions[0] / static_cast<float>(maxDim),
      volume_rendering_parameters.volume_dimensions[1] / static_cast<float>(maxDim),
      volume_rendering_parameters.volume_dimensions[2] / static_cast<float>(maxDim),
  };

  upload_volume_data(dimensions, data);
  upload_transferfunction_data(-0.1f, 1.2f);
}


// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::draw_volume(igl::opengl::ViewerCore& core) {
  //
  //  Setup
  //
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

//  glClearColor(0.f, 0.0f, 0.f, 1.0f);
//  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


  //
  //  Pre-Rendering
  //
  glBindVertexArray(bounding_box.vao);
  glUseProgram(bounding_box.program);

  Eigen::Matrix4f scaling = Eigen::Matrix4f::Zero();
  scaling.diagonal() << volume_rendering_parameters.normalized_volume_dimensions[0],
                        volume_rendering_parameters.normalized_volume_dimensions[1],
                        volume_rendering_parameters.normalized_volume_dimensions[2],
                        1.f;
  Eigen::Matrix4f model = core.model * scaling;

  glUniformMatrix4fv(BoundingBox::uniform_location.model_matrix, 1, GL_FALSE,
                     model.data());

  glUniformMatrix4fv(BoundingBox::uniform_location.view_matrix, 1, GL_FALSE,
                     core.view.data());

  glUniformMatrix4fv(BoundingBox::uniform_location.projection_matrix, 1, GL_FALSE,
                     core.proj.data());

  // Render entry points of bounding box
  glBindFramebuffer(GL_FRAMEBUFFER, bounding_box.entry_framebuffer);
  glClear(GL_COLOR_BUFFER_BIT);
  glCullFace(GL_FRONT);
  glDrawElements(GL_TRIANGLES, 12 * 3, GL_UNSIGNED_BYTE, nullptr);

  // Render exit points of bounding box
  glBindFramebuffer(GL_FRAMEBUFFER, bounding_box.exit_framebuffer);
  glClear(GL_COLOR_BUFFER_BIT);
  glCullFace(GL_BACK);
  glDrawElements(GL_TRIANGLES, 12 * 3, GL_UNSIGNED_BYTE, nullptr);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);


  //
  //  Volume rendering
  //
  glUseProgram(VolumeRendering::program);

  // Entry points texture
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, bounding_box.entry_texture);
  glUniform1i(VolumeRendering::uniform_location.entry_texture, 0);

  // Exit points texture
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, bounding_box.exit_texture);
  glUniform1i(VolumeRendering::uniform_location.entry_texture, 1);

  // Volume texture
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_3D, volume_rendering.volume_texture);
  glUniform1i(VolumeRendering::uniform_location.volume_texture, 2);

  // Transfer function texture
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_1D, volume_rendering.transfer_function_texture);
  glUniform1i(VolumeRendering::uniform_location.transfer_function, 3);

  glUniform1f(VolumeRendering::uniform_location.sampling_rate,
              volume_rendering_parameters.sampling_rate);


  // Rendering parameters
  glUniform3uiv(VolumeRendering::uniform_location.volume_dimensions, 3,
                volume_rendering_parameters.volume_dimensions.data());
  glUniform3fv(VolumeRendering::uniform_location.volume_dimensions_rcp, 3,
               volume_rendering_parameters.volume_dimensions_rcp.data());
  glUniform3fv(VolumeRendering::uniform_location.light_position, 3,
               core.light_position.data());
  glUniform3f(VolumeRendering::uniform_location.light_color_ambient, 0.5f, 0.5f, 0.5f);
  glUniform3f(VolumeRendering::uniform_location.light_color_diffuse, 0.8f, 0.8f, 0.8f);
  glUniform3f(VolumeRendering::uniform_location.light_color_specular, 1.f, 1.f, 1.f);
  glUniform1f(VolumeRendering::uniform_location.light_exponent_specular, 10.f);


  glDisable(GL_CULL_FACE);
  glDrawArrays(GL_TRIANGLES, 0, 6);

  glUseProgram(0);
  glBindVertexArray(0);

  return true;
}

// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::upload_volume_data(const Eigen::RowVector3i& tex_size,
                                               const Eigen::VectorXd& texture)
{
  std::vector<uint8_t> volume_data(texture.size());
  std::transform(
        texture.data(),
        texture.data() + texture.size(),
        volume_data.begin(),
        [](double d) {
    return static_cast<uint8_t>(d * std::numeric_limits<uint8_t>::max());
  }
  );

  glBindTexture(GL_TEXTURE_3D, volume_rendering.volume_texture);
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, tex_size[0], tex_size[1], tex_size[2], 0,
      GL_RED, GL_UNSIGNED_BYTE, volume_data.data());
  glBindTexture(GL_TEXTURE_3D, 0);
}


// TODO: Real transfer function
// TODO: Check opengl status codes and return false on failure
bool igl::opengl::VolumeGL::upload_transferfunction_data(float offset, float incline) {
  constexpr const int TransferFunctionWidth = 512;
  std::vector<std::array<uint8_t, 4>> transfer_function_data(512);

  // Create greyscale ramp
  for (int i = 0; i < TransferFunctionWidth; ++i) {
      float v = static_cast<float>(i * incline) / (TransferFunctionWidth - 1) + offset;
      if (v > 1.f) {
          v = 1.f;
      }
      if (v < 0.f) {
          v = 0.f;
      }
      const int val = v * std::numeric_limits<uint8_t>::max();
      const uint8_t value = static_cast<uint8_t>(val);

      transfer_function_data[i] = { value, value, value, value };
  }

  glBindTexture(GL_TEXTURE_1D, volume_rendering.transfer_function_texture);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, TransferFunctionWidth, 0, GL_RGBA,
               GL_UNSIGNED_BYTE, transfer_function_data.data());
  glBindTexture(GL_TEXTURE_1D, 0);
}


bool igl::opengl::VolumeGL::init_shaders() {

  // Shader transforming the vertices from model coordinates to clip space
  constexpr const char* BoundingBoxVertexShader = R"(
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
    )";

  // Using Krueger-Westermann rendering encodes the position of the vertex as its color
  constexpr const char* EntryBoxFragmentShader = R"(
    #version 150
      in vec3 color;
      out vec4 out_color;

      void main() {
        out_color = vec4(color, 1.0);
      }
    )";

  BoundingBox::program = igl::opengl::create_shader_program(
        BoundingBoxVertexShader,
        EntryBoxFragmentShader,
        { { "in_position", 0 } }
    );
  BoundingBox::uniform_location.model_matrix = glGetUniformLocation(
        BoundingBox::program, "model_matrix");
  BoundingBox::uniform_location.view_matrix = glGetUniformLocation(
        BoundingBox::program, "view_matrix");
  BoundingBox::uniform_location.projection_matrix = glGetUniformLocation(
        BoundingBox::program, "projection_matrix");

  cout << "BoundingBox program: " << BoundingBox::program << endl;
  cout << "BoundingBox uniform model_matrix: " << BoundingBox::uniform_location.model_matrix << endl;
  cout << "BoundingBox uniform view_matrix: " << BoundingBox::uniform_location.view_matrix << endl;
  cout << "BoundingBox uniform projection_matrix: " << BoundingBox::uniform_location.projection_matrix << endl;

  // Vertex shader that is used to trigger the volume rendering by rendering a static
  // screen-space filling quad.
  constexpr const char* VolumeRenderingVertexShader = R"(
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
      // Clipspace \in [-1, 1]
      gl_Position = vec4(positions[gl_VertexID], 0.0, 1.0);

      // uv coordinate s\in [0, 1]
      uv = (positions[gl_VertexID] + 1.0) / 2.0;
    })";

  // Shader that performs the actual volume rendering
  // Steps:
  // 1. Compute the ray direction by exit point color - entry point color
  // 2. Sample the volume along the ray
  // 3. Convert sample to color using the transfer function
  // 4. Compute central difference gradient
  // 5. Use the gradient for Phong shading
  // 6. Perform front-to-back compositing
  // 7. Stop if either the ray is exhausted or the combined transparency is above an
  //    early-ray termination threshold (0.99 in this case)
  constexpr const char* VolumeRenderingFragmentShader = R"(
    #version 150
    in vec2 uv;
    out vec4 out_color;

    uniform sampler2D entry_texture;
    uniform sampler2D exit_texture;

    uniform sampler3D volume_texture;
    uniform sampler1D transfer_function;

    uniform ivec3 volume_dimensions;
    uniform vec3 volume_dimensions_rcp;
    uniform float sampling_rate;

    struct Light_Parameters {
      vec3 position;
      vec3 ambient_color;
      vec3 diffuse_color;
      vec3 specular_color;
      float specular_exponent;
    };
    uniform Light_Parameters light_parameters;


    // Early-ray termination
    const float ERT_THRESHOLD = 0.99;
    const float REF_SAMPLING_INTERVAL = 150.0;

    vec3 centralDifferenceGradient(vec3 pos) {
      vec3 f;
      f.x = texture(volume_texture, pos + vec3(volume_dimensions_rcp.x, 0.0, 0.0)).a;
      f.y = texture(volume_texture, pos + vec3(0.0, volume_dimensions_rcp.y, 0.0)).a;
      f.z = texture(volume_texture, pos + vec3(0.0, 0.0, volume_dimensions_rcp.z)).a;

      vec3 b;
      b.x = texture(volume_texture, pos - vec3(volume_dimensions_rcp.x, 0.0, 0.0)).a;
      b.y = texture(volume_texture, pos - vec3(0.0, volume_dimensions_rcp.y, 0.0)).a;
      b.z = texture(volume_texture, pos - vec3(0.0, 0.0, volume_dimensions_rcp.z)).a;

      return (f - b) / 2.0;
    }

    vec3 blinn_phong(Light_Parameters light, vec3 material_ambient_color,
                     vec3 material_diffuse_color, vec3 material_specular_color,
                     vec3 position, vec3 normal, vec3 direction_to_camera)
    {
      vec3 direction_to_light = normalize(light.position - position);
      vec3 ambient = material_ambient_color * light.ambient_color;
      vec3 diffuse = material_diffuse_color * light.diffuse_color *
                     max(dot(normal, direction_to_light), 0.0);
      vec3 specular;
      {
        vec3 half_way_vector = normalize(direction_to_camera + direction_to_light);
        specular = material_specular_color * light.specular_color *
                   pow(max(dot(normal, half_way_vector), 0.0), light.specular_exponent);
      }

      return ambient + diffuse + specular;
    }

    void main() {
      vec3 entry = texture(entry_texture, uv).rgb;
      vec3 exit = texture(exit_texture, uv).rgb;
      if (entry == exit) {
        out_color = vec4(0.0, 0.0, 0.0, 0.0);
        return;
      }

      // Combined final color that the volume rendering computed
      vec4 result = vec4(0.0);

      vec3 ray_direction = exit - entry;

      float t_end = length(ray_direction);
      float t_incr = min(
        t_end,
        t_end / (sampling_rate * length(ray_direction * volume_dimensions))
      );
      t_incr = 0.01;

      vec3 normalized_ray_direction = normalize(ray_direction);

      float t = 0.0;
      while (t < t_end) {
        vec3 sample_pos = entry + t * normalized_ray_direction;
        float value = texture(volume_texture, sample_pos).r;
        vec4 color = texture(transfer_function, value);
        if (color.a > 0) {
          // Gradient
          vec3 gradient = centralDifferenceGradient(sample_pos);

          // Lighting
          //color.rgb = blinn_phong(light_parameters, color.rgb, color.rgb, vec3(1.0),
                                  //sample_pos, gradient, -normalized_ray_direction);

          // Front-to-back Compositing
          color.a = 1.0 - pow(1.0 - color.a, t_incr * REF_SAMPLING_INTERVAL);
          result.rgb = result.rgb + (1.0 - result.a) * color.a * color.rgb;
          result.a = result.a + (1.0 - result.a) * color.a;
        }

        if (result.a > ERT_THRESHOLD) {
          t = t_end;
        }
        else {
          t += t_incr;
        }
      }

      result.a = length(result.rgb);
      out_color = result;
    })";
  VolumeRendering::program = igl::opengl::create_shader_program(
        VolumeRenderingVertexShader,
        VolumeRenderingFragmentShader, {});

  VolumeRendering::uniform_location.entry_texture = glGetUniformLocation(
      VolumeRendering::program, "entry_texture");
  VolumeRendering::uniform_location.exit_texture = glGetUniformLocation(
      VolumeRendering::program, "exit_texture");
  VolumeRendering::uniform_location.volume_texture = glGetUniformLocation(
      VolumeRendering::program, "volume_texture");
  VolumeRendering::uniform_location.volume_dimensions = glGetUniformLocation(
      VolumeRendering::program, "volume_dimensions");
  VolumeRendering::uniform_location.volume_dimensions_rcp = glGetUniformLocation(
      VolumeRendering::program, "volume_dimensions_rcp");
  VolumeRendering::uniform_location.transfer_function = glGetUniformLocation(
      VolumeRendering::program, "transfer_function");
  VolumeRendering::uniform_location.sampling_rate = glGetUniformLocation(
      VolumeRendering::program, "sampling_rate");
  VolumeRendering::uniform_location.light_position = glGetUniformLocation(
      VolumeRendering::program, "light_parameters.position");
  VolumeRendering::uniform_location.light_color_ambient = glGetUniformLocation(
      VolumeRendering::program, "light_parameters.ambient_color");
  VolumeRendering::uniform_location.light_color_diffuse = glGetUniformLocation(
      VolumeRendering::program, "light_parameters.diffuse_color");
  VolumeRendering::uniform_location.light_color_specular = glGetUniformLocation(
      VolumeRendering::program, "light_parameters.specular_color");
  VolumeRendering::uniform_location.light_exponent_specular = glGetUniformLocation(
      VolumeRendering::program, "light_parameters.specular_exponent");


  cout << "VolumeRendering program: " << VolumeRendering::program << endl;
  cout << "VolumeRendering uniform entry_texture: " << VolumeRendering::uniform_location.entry_texture << endl;
  cout << "VolumeRendering uniform exit_texture: " << VolumeRendering::uniform_location.exit_texture << endl;
  cout << "VolumeRendering uniform volume_texture: " << VolumeRendering::uniform_location.volume_texture << endl;
  cout << "VolumeRendering uniform volume_dimensions: " << VolumeRendering::uniform_location.volume_dimensions << endl;
  cout << "VolumeRendering uniform volume_dimensions_rcp: " << VolumeRendering::uniform_location.volume_dimensions_rcp << endl;
  cout << "VolumeRendering uniform transfer_function: " << VolumeRendering::uniform_location.transfer_function << endl;
  cout << "VolumeRendering uniform sampling_rate: " << VolumeRendering::uniform_location.sampling_rate << endl;
  cout << "VolumeRendering uniform light_position: " << VolumeRendering::uniform_location.light_position << endl;
  cout << "VolumeRendering uniform light_color_ambient: " << VolumeRendering::uniform_location.light_color_ambient << endl;
  cout << "VolumeRendering uniform light_color_diffuse: " << VolumeRendering::uniform_location.light_color_diffuse << endl;
  cout << "VolumeRendering uniform light_color_specular: " << VolumeRendering::uniform_location.light_color_specular << endl;
  cout << "VolumeRendering uniform light_exponent_specular: " << VolumeRendering::uniform_location.light_exponent_specular << endl;
}
