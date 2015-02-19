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
#include "sys/intrinsics.h"
#include "sys/sysinfo.h"
#include "sys/sync/atomic.h"
#include "sys/stl/vector.h"
#include "sys/stl/string.h"

#include "math/math.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include "math/bbox.h"
#include "math/affinespace.h"

#include "simd/simd.h"

  /*! Ray structure. Contains all information about a ray including
   *  precomputed reciprocal direction. */
  struct RTCRay
  {
    /*! Default construction does nothing. */
    __forceinline RTCRay() {}

    /*! Constructs a ray from origin, direction, and ray segment. Near
     *  has to be smaller than far. */
    __forceinline RTCRay(const embree::Vec3fa& org, const embree::Vec3fa& dir, 
			 float tnear = embree::zero, float tfar = embree::inf, 
			 float time = embree::zero, int mask = -1)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), geomID(-1), primID(-1), instID(-1), mask(mask), time(time) {}

    /*! Tests if we hit something. */
    __forceinline operator bool() const { return geomID != -1; }

  public:
    embree::Vec3fa org;        //!< Ray origin
    embree::Vec3fa dir;        //!< Ray direction
    float tnear;       //!< Start of ray segment
    float tfar;        //!< End of ray segment
    float time;        //!< Time of this ray for motion blur.
    int mask;          //!< used to mask out objects during traversal

  public:
    embree::Vec3fa Ng;         //!< Not normalized geometry normal
    float u;           //!< Barycentric u coordinate of hit
    float v;           //!< Barycentric v coordinate of hit
    int geomID;           //!< geometry ID
    int primID;           //!< primitive ID
    int instID;           //!< instance ID
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const RTCRay& ray) {
    return cout << "{ " << 
      "org = " << ray.org << ", dir = " << ray.dir << ", near = " << ray.tnear << ", far = " << ray.tfar << ", time = " << ray.time << ", " <<
      "instID = " << ray.instID <<  ", geomID = " << ray.geomID << ", primID = " << ray.primID <<  ", " << "u = " << ray.u <<  ", v = " << ray.v << ", Ng = " << ray.Ng << " }";
}
