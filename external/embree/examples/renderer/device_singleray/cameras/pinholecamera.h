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

#ifndef __EMBREE_PINHOLE_CAMERA_H__
#define __EMBREE_PINHOLE_CAMERA_H__

#include "camera.h"

namespace embree
{
  /*! Implements the pinhole camera model. */
  class PinHoleCamera : public Camera
  {
  public:

    /*! Construction from parameter container. */
    PinHoleCamera(const Parms& parms) {
      local2world = parms.getTransform("local2world");
      angle       = parms.getFloat("angle",64.0f);
      aspectRatio = parms.getFloat("aspectRatio",1.0f);
      Vector3f W     = xfmVector(local2world, Vector3f(-0.5f*aspectRatio,-0.5f,0.5f*rcp(tanf(deg2rad(0.5f*angle)))));
      pixel2world = AffineSpace3f(aspectRatio*local2world.l.vx,local2world.l.vy,W,local2world.p);
    }

    void ray(const Vec2f& pixel, const Vec2f& sample, Ray& ray_o) const {
      new (&ray_o) Ray(pixel2world.p,normalize(pixel.x*pixel2world.l.vx + (1.0f-pixel.y)*pixel2world.l.vy + pixel2world.l.vz));
    }

  protected:
    float aspectRatio;
    AffineSpace3f local2world;    //!< transformation from camera space to world space
    AffineSpace3f pixel2world;    //!< special transformation to generate rays
  };
}

#endif
