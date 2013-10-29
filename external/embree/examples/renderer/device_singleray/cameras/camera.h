// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifndef __EMBREE_CAMERA_H__
#define __EMBREE_CAMERA_H__

#include "../api/parms.h"
#include "embree/common/ray.h"

namespace embree
{
  /*! Interface to different camera models. */
  class Camera : public RefCount {
    ALIGNED_CLASS
  public:

    /*! Virtual interface destructor. */
    virtual ~Camera() {}

    /*! Computes a primary ray for a pixel. */
    virtual void ray(const Vec2f& pixel,  /*!< The pixel location on the on image plane in the range from 0 to 1. */
                     const Vec2f& sample, /*!< The lens sample in [0,1) for depth of field. */
                     Ray& ray_o)          /*!< To return the ray. */ const = 0;

    /*! Field of view. */
    float angle;
  };
}

#endif
