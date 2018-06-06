// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Francis Williams <francis@fwilliams.info>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VOLUMEDATA_H
#define VOLUMEDATA_H

#include "../igl_inline.h"
#include "VolumeGL.h"

#include <Eigen/Core>


namespace igl {
namespace opengl {

class VolumeData
{
public:
  VolumeData(int id) : _id(id) {}

  IGL_INLINE void set_data(const Eigen::RowVector3i& size, const Eigen::VectorXd& voldata);

  // TODO
  IGL_INLINE void set_transfer_function();

  IGL_INLINE void draw(igl::opengl::ViewerCore& core);

  ~VolumeData() {
    volumegl.free();
  }

  int id() const {
    return _id;
  }

  float sampling_rate = 128.0;

private:
  enum DirtyFlags
  {
    DIRTY_NONE              = 0x0000,
    DIRTY_VOLUME            = 0x0001,
    DIRTY_TRANSFER_FUNCTION = 0x0002,
    DIRTY_ALL               = 0xFFFF
  };

  // What we think the viewport size is. If it changes, we re-allocate textures.
  Eigen::Vector4f viewport_size = {-1.0, -1.0, -1.0, -1.0};
  Eigen::RowVector3i size;
  Eigen::VectorXd voldata;

  igl::opengl::VolumeGL volumegl;

  // Marks dirty buffers that need to be uploaded to OpenGL
  uint32_t dirty;

  // Unique identifier
  int _id;
};

} // namespace opengl
} // namespace igl

#endif // VOLUMEDATA_H
