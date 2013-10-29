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

#ifndef __EMBREE_DEPTH_OF_FIELD_CAMERA_H__
#define __EMBREE_DEPTH_OF_FIELD_CAMERA_H__

#include "pinholecamera.h"

namespace embree
{
  /*! Implements a depth of field camera model. */
  class DepthOfFieldCamera : public PinHoleCamera
  {
  public:

    /*! Construction from parameter container. */
    DepthOfFieldCamera (const Parms& parms) : PinHoleCamera(parms) {
      lensRadius     = parms.getFloat("lensRadius",0.0f);
      focalDistance  = parms.getFloat("focalDistance");
      focalDistance /= length(0.5f*pixel2world.l.vx + 0.5f*pixel2world.l.vy + pixel2world.l.vz);
    }

    void ray(const Vec2f& pixel, const Vec2f& sample, Ray& ray_o) const 
    {
      Vec2f lens = uniformSampleDisk(sample,lensRadius);
      Vector3f begin = xfmPoint(local2world, Vector3f(lens.x,lens.y,0.0f));
      Vector3f end   = pixel2world.p + focalDistance*(pixel.x*pixel2world.l.vx + (1.0f-pixel.y)*pixel2world.l.vy + pixel2world.l.vz);
      new (&ray_o) Ray(begin, normalize(end - begin));
    }

  protected:
    float lensRadius;     //!< radius of the lens
    float focalDistance;  //!< distance of the focal plane
  };
}

#endif
