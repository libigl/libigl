// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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
  class PinholeCamera : public Camera
  {
  public:

    /*! Construction from parameter container. */
    PinholeCamera(const Parms& parms) {
      localToWorld      = parms.getTransform("local2world");
      angle             = parms.getFloat("angle",64.0f);
      float aspectRatio = parms.getFloat("aspectRatio",1.0f);
      Vec3f W           = xfmVector(localToWorld, Vec3f(-0.5f*aspectRatio,-0.5f,0.5f*rcp(tanf(deg2rad(0.5f*angle)))));
      xfm = AffineSpace3f(aspectRatio*localToWorld.l.vx,localToWorld.l.vy,W,localToWorld.p);
    }

    void ray(const Vec2f& pixel, const Vec2f& sample, Ray& ray_o) const {
      new (&ray_o) Ray(xfm.p,normalize(pixel.x*xfm.l.vx + (1.0f-pixel.y)*xfm.l.vy + xfm.l.vz));
    }

  protected:
    AffineSpace3f localToWorld; //!< transformation from camera space to world space
    AffineSpace3f xfm;          //!< special transformation to generate rays
  };
}

#endif
