#ifndef IGL_OPENGL_GFLW_IMGUI_SELECTIONWIDGET_H
#define IGL_OPENGL_GFLW_IMGUI_SELECTIONWIDGET_H
#include "../../../igl_inline.h"
#include "ImGuiWidget.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <functional>

namespace igl{ namespace opengl{ namespace glfw{ namespace imgui{

class SelectionWidget: public ImGuiWidget
{
public:
  // customizable hotkeys
  std::string MARQUEE_KEY = "Mm";
  // leave 'L' for show_lines in viewer
  std::string LASSO_KEY = "l";
  std::string OFF_KEY = "Vv";
  enum Mode
  {
    OFF                 = 0,
    RECTANGULAR_MARQUEE = 1,
    ELLIPTICAL_MARQUEE  = 2,
    POLYGONAL_LASSO     = 3,
    LASSO               = 4,
    NUM_MODES           = 5
  } mode = RECTANGULAR_MARQUEE;
  bool is_down = false;
  bool has_moved_since_down = false;
  bool is_drawing = false;
  // min and max corners of 2D rectangular marquee
  Eigen::Matrix<float,2,2> M = Eigen::Matrix<float,2,2>::Zero();
  // list of points of 2D lasso marquee
  std::vector<Eigen::RowVector2f> L;
  // callback called when slection is completed (usually on mouse_up)
  std::function<void(void)> callback;
  // callback called after mode is changed
  std::function<void(Mode)> callback_post_mode_change;
  IGL_INLINE virtual void init(Viewer *_viewer, ImGuiPlugin *_plugin) override;
  IGL_INLINE virtual void draw() override;
  IGL_INLINE virtual bool mouse_down(int button, int modifier) override;
  IGL_INLINE virtual bool mouse_up(int button, int modifier) override;
  IGL_INLINE virtual bool mouse_move(int mouse_x, int mouse_y) override;
  IGL_INLINE virtual bool key_pressed(unsigned int key, int modifiers) override;
  IGL_INLINE void clear();
  // helpers
  IGL_INLINE static void circle(const Eigen::Matrix<float,2,2> & M,  std::vector<Eigen::RowVector2f> & L);
  IGL_INLINE static void rect(const Eigen::Matrix<float,2,2> & M,  std::vector<Eigen::RowVector2f> & L);
  IGL_INLINE static Eigen::RowVector2f xy(const Viewer * v);
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}}}}

#ifndef IGL_STATIC_LIBRARY
#include "SelectionWidget.cpp"
#endif
#endif
