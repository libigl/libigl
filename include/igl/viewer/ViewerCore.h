#ifndef IGL_VIEWER_CORE_H
#define IGL_VIEWER_CORE_H

#include <igl/igl_inline.h>
#include <igl/viewer/TextRenderer.h>

namespace igl
{

class ViewerCore
#ifdef ENABLE_XML_SERIALIZATION
: public ::igl::XMLSerialization
#endif
{
public:
  IGL_INLINE ViewerCore();

  IGL_INLINE void init();
  IGL_INLINE void shut();


  IGL_INLINE void InitSerialization();

  IGL_INLINE void align_camera_center(const Eigen::MatrixXd& V); // Adjust the view to see the entire model

  // Determines how much to zoom and shift such that the mesh fills the unit
  // box (centered at the origin)
  IGL_INLINE void get_scale_and_shift_to_fit_mesh(
    const Eigen::MatrixXd& V,
    float & zoom,
    Eigen::Vector3f& shift);

  IGL_INLINE void clear_framebuffers();

  // Draw everything
  IGL_INLINE void draw(ViewerData& data, OpenGL_state& opengl);

  TextRenderer textrenderer;

  // Shape material
  float shininess;

  // Colors
  Eigen::Vector3f background_color;
  Eigen::Vector3f line_color;

  // Lighting
  Eigen::Vector3f light_position;
  float lighting_factor;

  // Trackball angle (quaternion)
  Eigen::Vector4f trackball_angle;

  // Model viewing parameters
  float model_zoom;
  Eigen::Vector3f model_translation;

  // Model viewing paramters (uv coordinates)
  float model_zoom_uv;
  Eigen::Vector3f model_translation_uv;

  // Camera parameters
  float camera_zoom;
  bool orthographic;
  Eigen::Vector3f camera_eye;
  Eigen::Vector3f camera_up;
  Eigen::Vector3f camera_center;
  float camera_view_angle;
  float camera_dnear;
  float camera_dfar;

  // Visualization options
  bool show_overlay;
  bool show_overlay_depth;
  bool show_texture;
  bool show_faces;
  bool show_lines;
  bool show_vertid;
  bool show_faceid;
  bool invert_normals;

  // Point size / line width
  float point_size;
  float line_width;

  // Animation
  bool is_animating;
  double animation_max_fps;

  // Caches the two-norm between the min/max point of the bounding box
  float object_scale;

  // Window size
  Eigen::Vector4f viewport;

  // Save the OpenGL transformation matrices used for the previous rendering pass
  Eigen::Matrix4f view;
  Eigen::Matrix4f model;
  Eigen::Matrix4f proj;
};

}

#ifndef IGL_STATIC_LIBRARY
#  include "ViewerCore.cpp"
#endif

#endif
