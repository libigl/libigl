#ifndef IGL_OPENGL_STATE_H
#define IGL_OPENGL_STATE_H

// Coverts mesh data inside a igl::ViewerData class in an OpenGL compatible format
// The class includes a shader and the opengl calls to plot the data

#include <igl/igl_inline.h>
#include <igl/viewer/OpenGL_shader.h>

namespace igl
{

class OpenGL_state
{
public:
  typedef unsigned int GLuint;

  GLuint vao_mesh;
  GLuint vao_overlay_lines;
  GLuint vao_overlay_points;
  OpenGL_shader shader_mesh;
  OpenGL_shader shader_overlay_lines;
  OpenGL_shader shader_overlay_points;

  GLuint vbo_V; // Vertices of the current mesh (#V x 3)
  GLuint vbo_V_uv; // UV coordinates for the current mesh (#V x 2)
  GLuint vbo_V_normals; // Vertices of the current mesh (#V x 3)
  GLuint vbo_V_ambient; // Ambient material  (#V x 3)
  GLuint vbo_V_diffuse; // Diffuse material  (#V x 3)
  GLuint vbo_V_specular; // Specular material  (#V x 3)

  GLuint vbo_F; // Faces of the mesh (#F x 3)
  GLuint vbo_tex; // Texture

  GLuint vbo_lines_F;         // Indices of the line overlay
  GLuint vbo_lines_V;         // Vertices of the line overlay
  GLuint vbo_lines_V_colors;  // Color values of the line overlay
  GLuint vbo_points_F;        // Indices of the point overlay
  GLuint vbo_points_V;        // Vertices of the point overlay
  GLuint vbo_points_V_colors; // Color values of the point overlay

  // Temporary copy of the content of each VBO
  Eigen::MatrixXf V_vbo;
  Eigen::MatrixXf V_normals_vbo;
  Eigen::MatrixXf V_ambient_vbo;
  Eigen::MatrixXf V_diffuse_vbo;
  Eigen::MatrixXf V_specular_vbo;
  Eigen::MatrixXf V_uv_vbo;
  Eigen::MatrixXf lines_V_vbo;
  Eigen::MatrixXf lines_V_colors_vbo;
  Eigen::MatrixXf points_V_vbo;
  Eigen::MatrixXf points_V_colors_vbo;

  int tex_u;
  int tex_v;
  Eigen::Matrix<char,Eigen::Dynamic,1> tex;

  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> F_vbo;
  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> lines_F_vbo;
  Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> points_F_vbo;

  // Marks dirty buffers that need to be uploaded to OpenGL
  uint32_t dirty;

  // Initialize shaders and buffers
  IGL_INLINE void init();

  // Release all resources
  IGL_INLINE void free();

  // Create a new set of OpenGL buffer objects
  IGL_INLINE void init_buffers();

  // Update contents from a 'Data' instance
  IGL_INLINE void set_data(const igl::ViewerData &data, bool face_based, bool invert_normals);

  // Bind the underlying OpenGL buffer objects for subsequent mesh draw calls
  IGL_INLINE void bind_mesh();

  /// Draw the currently buffered mesh (either solid or wireframe)
  IGL_INLINE void draw_mesh(bool solid);

  // Bind the underlying OpenGL buffer objects for subsequent line overlay draw calls
  IGL_INLINE void bind_overlay_lines();

  /// Draw the currently buffered line overlay
  IGL_INLINE void draw_overlay_lines();

  // Bind the underlying OpenGL buffer objects for subsequent point overlay draw calls
  IGL_INLINE void bind_overlay_points();

  /// Draw the currently buffered point overlay
  IGL_INLINE void draw_overlay_points();

  // Release the OpenGL buffer objects
  IGL_INLINE void free_buffers();

  enum DirtyFlags
  {
    DIRTY_NONE           = 0x0000,
    DIRTY_POSITION       = 0x0001,
    DIRTY_UV             = 0x0002,
    DIRTY_NORMAL         = 0x0004,
    DIRTY_AMBIENT        = 0x0008,
    DIRTY_DIFFUSE        = 0x0010,
    DIRTY_SPECULAR       = 0x0020,
    DIRTY_TEXTURE        = 0x0040,
    DIRTY_FACE           = 0x0080,
    DIRTY_MESH           = 0x00FF,
    DIRTY_OVERLAY_LINES  = 0x0100,
    DIRTY_OVERLAY_POINTS = 0x0200,
    DIRTY_ALL            = 0x03FF
  };

};

}

#ifndef IGL_STATIC_LIBRARY
#  include "OpenGL_state.cpp"
#endif

#endif
