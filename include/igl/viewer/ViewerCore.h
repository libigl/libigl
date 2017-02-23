// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_VIEWER_VIEWER_CORE_H
#define IGL_VIEWER_VIEWER_CORE_H

#ifdef IGL_VIEWER_WITH_NANOGUI
#include <igl/viewer/TextRenderer.h>
#endif
#include <igl/viewer/ViewerData.h>
#include <igl/viewer/OpenGL_state.h>

#include <igl/igl_inline.h>
#include <Eigen/Geometry>
#include <Eigen/Core>

namespace igl
{
namespace viewer
{

// Basic class of the 3D mesh viewer
// TODO: write documentation

class ViewerCore
{
public:
  IGL_INLINE ViewerCore();

  // Initialization
  IGL_INLINE void init();

  // Shutdown
  IGL_INLINE void shut();

  // ------------------- Camera control functions

  // Set camera center to new position (keep orientation)
  IGL_INLINE void set_camera_position(
    const Eigen::Vector3f& pos);

  // Adjust the camera to see the entire model
  IGL_INLINE void align_camera_center(
    const ViewerData& data);
  IGL_INLINE void align_camera_center(
    const Eigen::MatrixXd& V);
  IGL_INLINE void align_camera_center(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F);

  // Determines how much to zoom and shift such that the mesh fills the unit
  // box (centered at the origin)
  IGL_INLINE void get_zoom_and_shift_to_fit_mesh(
    const Eigen::MatrixXd& V,
    float & zoom,
    Eigen::Vector3f& shift);
  IGL_INLINE void get_zoom_and_shift_to_fit_mesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    float & zoom,
    Eigen::Vector3f& shift);

  // ------------------- Drawing functions

  // Clear the frame buffers
  IGL_INLINE void clear_framebuffers();

  // Draw everything
  IGL_INLINE void draw(ViewerData& data, OpenGL_state& opengl, bool update_matrices = true);
  
  IGL_INLINE void draw_buffer(
    ViewerData& data,
    OpenGL_state& opengl,
    bool update_matrices,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A);

  IGL_INLINE void draw_buffer(
    std::vector<ViewerData*>& data,
    std::vector<OpenGL_state*>& opengl,
    bool update_matrices,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A);

  // Trackball angle (quaternion)
  enum RotationType
  {
    ROTATION_TYPE_TRACKBALL = 0,
    ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
    NUM_ROTATION_TYPES = 2
  };
  IGL_INLINE void set_rotation_type(const RotationType & value);

  // ------------------- Properties

#ifdef IGL_VIEWER_WITH_NANOGUI
  // Text rendering helper
  TextRenderer textrenderer;
#endif

  // Colors
  Eigen::Vector4f background_color;
  Eigen::Vector4f line_color;

  // Lighting
  float shininess;
  Eigen::Vector3f light_position;
  float lighting_factor;

  // Global scene transformation
  RotationType rotation_type;
  Eigen::Quaternionf trackball_angle;
  Eigen::Vector3f global_translation;

  // Camera parameters
  float camera_zoom;
  bool orthographic;
  Eigen::Vector3f camera_eye;
  Eigen::Vector3f camera_up;
  Eigen::Vector3f camera_center;
  float camera_view_angle;
  float camera_dnear;
  float camera_dfar;

  // Animation
  bool is_animating;
  double animation_max_fps;

  // Viewport size
  Eigen::Vector4f viewport;

  // Save the OpenGL transformation matrices used for the previous rendering pass
  Eigen::Matrix4f view;
  Eigen::Matrix4f proj;
  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}
}

#ifndef IGL_STATIC_LIBRARY
#  include "ViewerCore.cpp"
#endif

#endif
