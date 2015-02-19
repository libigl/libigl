// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "sys/platform.h"
#include "sys/ref.h"
#include "math/math.h"
#include "math/vec3.h"
//#include "math/color.h"
#include "math/affinespace.h"

namespace embree
{
  /* camera settings */
  struct Camera 
  {
  public:
    Camera () 
      : from(0.0001f,0.0001f,-3.0f), to(0,0,0), up(0,1,0), fov(90) {}

    Camera (Vec3fa& from, Vec3fa& to, Vec3fa& up, float fov) 
      : from(from), to(to), up(up), fov(fov) {}

    AffineSpace3fa camera2world () { return AffineSpace3fa::lookat(from, to, up); }
    AffineSpace3fa world2camera () { return rcp(AffineSpace3fa::lookat(from, to, up)); }
    Vec3fa world2camera(const Vec3fa& p) { return xfmPoint(world2camera(),p); }
    Vec3fa camera2world(const Vec3fa& p) { return xfmPoint(camera2world(),p); }

    AffineSpace3fa pixel2world (size_t width, size_t height) 
    {
      const float fovScale = 1.0f/tanf(deg2rad(0.5f*fov));
      const AffineSpace3fa local2world = AffineSpace3fa::lookat(from, to, up);
      return AffineSpace3fa(local2world.l.vx,
                           -local2world.l.vy,
                           -0.5f*width*local2world.l.vx + 0.5f*height*local2world.l.vy + 0.5f*height*fovScale*local2world.l.vz,
                           local2world.p);
    }

    void move (float dx, float dy, float dz)
    {
      AffineSpace3fa xfm = camera2world();
      Vec3fa ds = xfmVector(xfm,Vec3fa(dx,dy,dz));
      from += ds;
      to   += ds;
    }

    void rotate (float dtheta, float dphi)
    {
      Vec3fa view = xfmPoint(world2camera(),to);
      float theta = atan2f(view.x, view.z); theta += dtheta;
      float phi   = asinf (view.y);         phi   += dphi;
      float x = cosf(phi)*sinf(theta);
      float y = sinf(phi);
      float z = cosf(phi)*cosf(theta);
      to = xfmPoint(camera2world(),length(view)*Vec3fa(x,y,z));
    }

    void rotateOrbit (float dtheta, float dphi)
    {
      Vec3fa view = normalize(xfmVector(world2camera(),to - from));
      float theta = atan2f(view.x, view.z); theta += dtheta;
      float phi   = asinf (view.y);         phi   += dphi;
      float x = cosf(phi)*sinf(theta);
      float y = sinf(phi);
      float z = cosf(phi)*cosf(theta);
      Vec3fa view1 = xfmVector(camera2world(),Vec3fa(x,y,z));
      from = to - length(to - from) * view1;
    }

    void dolly (float ds)
    {
      float dollySpeed = 0.01f;
      float k = powf((1.0f-dollySpeed), ds);
      from += length(to-from) * (1-k) * normalize(to-from);
    }

  public:
    Vec3fa from;   //!< position of camera
    Vec3fa to;     //!< look at point
    Vec3fa up;     //!< up vector
    float fov;       //!< field of view
  };
}
