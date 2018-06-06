#ifndef VOLUMEGL_H
#define VOLUMEGL_H

#include "gl.h"
#include "create_shader_program.h"

#include <Eigen/Core>
#include <array>

namespace igl {
namespace opengl {

class ViewerCore;

class VolumeGL
{
public:
  // TODO: Check opengl status codes and return false on failure
  bool init(igl::opengl::ViewerCore& core);

  void free();

  // TODO: Check opengl status codes and return false on failure
  bool set_data(const Eigen::RowVector3i& dimensions, const Eigen::VectorXd& data);

  bool draw(igl::opengl::ViewerCore& core, GLfloat sampling_rate);

  // TODO: Check opengl status codes and return false on failure
  bool resize_framebuffer_textures(ViewerCore& core);

  bool is_initialized() const
  {
    return _is_initialized;
  }

  Eigen::Vector3i volume_dimensions() const
  {
    return Eigen::Vector3i(
        volume_rendering_parameters.volume_dimensions[0],
        volume_rendering_parameters.volume_dimensions[1],
        volume_rendering_parameters.volume_dimensions[2]);
  }

  // TODO: Check opengl status codes and return false on failure
  static bool init_shaders();

private:
  // TODO: Check opengl status codes and return false on failure
  bool upload_volume_data(const Eigen::RowVector3i& tex_size,
                          const Eigen::VectorXd& texture);

  // TODO: Real transfer function
  // TODO: Check opengl status codes and return false on failure
  bool upload_transferfunction_data(float offset, float incline);

  bool _is_initialized = false;

  struct BoundingBox
  {
    // The shader program gets initialized once and for all when init_shaders() gets called
    static GLuint program;

    GLuint vao;
    GLuint vbo;
    GLuint ibo;

    GLuint entry_framebuffer;
    GLuint entry_texture;

    GLuint exit_framebuffer;
    GLuint exit_texture;

    struct UniformLocation {
      GLint model_matrix = -1;
      GLint view_matrix = -1;
      GLint projection_matrix = -1;
    };
    static UniformLocation uniform_location;
  } bounding_box;

  struct VolumeRendering
  {
    // The shader program gets initialized once and for all when init_shaders() gets called
    static GLuint program;

    GLuint volume_texture;
    GLuint transfer_function_texture;

    struct Uniform_Location {
      GLint entry_texture = -1;
      GLint exit_texture = -1;
      GLint volume_texture = -1;
      GLint transfer_function = -1;

      GLint volume_dimensions = -1;
      GLint volume_dimensions_rcp = -1;
      GLint sampling_rate = -1;

      GLint light_position = -1;
      GLint light_color_ambient = -1;
      GLint light_color_diffuse = -1;
      GLint light_color_specular = -1;
      GLint light_exponent_specular = -1;
    };
    static Uniform_Location uniform_location;
  } volume_rendering;

  struct VolumeRenderingParameters
  {
    std::array<GLuint, 3> volume_dimensions;
    std::array<GLfloat, 3> volume_dimensions_rcp;
    std::array<float, 3> normalized_volume_dimensions;
  } volume_rendering_parameters;
};

} // namespace opengl
} // namespace igl

#endif // VOLUMEGL_H
