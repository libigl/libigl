#include "VolumeData.h"
#include "ViewerCore.h"

#include <iostream>

using namespace std;


IGL_INLINE void igl::opengl::VolumeData::set_data(const Eigen::RowVector3i& size, const Eigen::VectorXd& voldata) {
  this->size = size;
  this->voldata = voldata;
  dirty |= VolumeData::DIRTY_VOLUME;
}

IGL_INLINE void igl::opengl::VolumeData::set_transfer_function() {}


IGL_INLINE void igl::opengl::VolumeData::draw(igl::opengl::ViewerCore& core) {
  if (!volumegl.is_initialized()) {
    cout << "Volume not initialized, initializing volume" << endl;
    volumegl.init(core);
  }

  if (dirty & VolumeData::DIRTY_VOLUME) {
    cout << "Volume data is dirty. Setting it." << endl;
    volumegl.set_data(size, voldata);
  }

  if (dirty & VolumeData::DIRTY_TRANSFER_FUNCTION) {
    // TODO: Set the transfer function
  }

  if (core.viewport[0] != viewport_size[0] || core.viewport[1] != viewport_size[1] || core.viewport[2] != viewport_size[2] || core.viewport[3] != viewport_size[3]) {
    volumegl.resize_framebuffer_textures(core.viewport);
    viewport_size = core.viewport;
  }

  dirty = VolumeData::DIRTY_NONE;

  volumegl.draw_volume(core);
}
