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

#ifndef __EMBREE_TUTORIALS_CAMERA_H__
#define __EMBREE_TUTORIALS_CAMERA_H__

#include "sys/platform.h"
#include "sys/ref.h"
#include "math/math.h"
#include "math/vec3.h"
#include "math/color.h"
#include "math/affinespace.h"

namespace embree
{
  /* camera settings */
  struct Camera 
  {
  public:
    Camera () 
      : from(0.0001,0.0001,-3), to(0,0,0), up(0,1,0), fov(90) {}

    Camera (Vector3f from, Vector3f to, Vector3f up, float fov) 
      : from(from), to(to), up(up), fov(fov) {}

    AffineSpace3f camera2world () { return AffineSpace3f::lookAtPoint(from, to, up); }
    AffineSpace3f world2camera () { return rcp(AffineSpace3f::lookAtPoint(from, to, up)); }
    Vector3f world2camera(const Vector3f& p) { return xfmPoint(world2camera(),p); }
    Vector3f camera2world(const Vector3f& p) { return xfmPoint(camera2world(),p); }

    AffineSpace3f pixel2world (size_t width, size_t height) 
    {
      const float fovScale = 1.0f/tanf(deg2rad(0.5f*fov));
      const AffineSpace3f local2world = AffineSpace3f::lookAtPoint(from, to, up);
      return AffineSpace3f(local2world.l.vx,
                           -local2world.l.vy,
                           -0.5f*width*local2world.l.vx + 0.5f*height*local2world.l.vy + 0.5f*height*fovScale*local2world.l.vz,
                           local2world.p);
    }

    void move (float dx, float dy, float dz)
    {
      AffineSpace3f xfm = camera2world();
      Vector3f ds = xfmVector(xfm,Vector3f(dx,dy,dz));
      from += ds;
      to   += ds;
    }

    void rotate (float dtheta, float dphi)
    {
      Vector3f view = xfmPoint(world2camera(),to);
      float theta = atan2f(view.x, view.z); theta += dtheta;
      float phi   = asinf (view.y);         phi   += dphi;
      float x = cosf(phi)*sinf(theta);
      float y = sinf(phi);
      float z = cosf(phi)*cosf(theta);
      to = xfmPoint(camera2world(),length(view)*Vector3f(x,y,z));
    }

    void rotateOrbit (float dtheta, float dphi)
    {
      Vector3f view = normalize(xfmVector(world2camera(),to - from));
      float theta = atan2f(view.x, view.z); theta += dtheta;
      float phi   = asinf (view.y);         phi   += dphi;
      float x = cosf(phi)*sinf(theta);
      float y = sinf(phi);
      float z = cosf(phi)*cosf(theta);
      Vector3f view1 = xfmVector(camera2world(),Vector3f(x,y,z));
      from = to - length(to - from) * view1;
    }

    void dolly (float ds)
    {
      float dollySpeed = 0.01f;
      float k = powf((1.0f-dollySpeed), ds);
      from += length(to-from) * (1-k) * normalize(to-from);
    }

  public:
    Vector3f from;   //!< position of camera
    Vector3f to;     //!< look at point
    Vector3f up;     //!< up vector
    float fov;       //!< field of view
  };
}

#endif
