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

#ifndef __EMBREE_RAY4_H__
#define __EMBREE_RAY4_H__

#include "ray.h"

namespace embree
{
  /*! Ray structure. Contains all information of 4 rays including
   *  precomputed reciprocal direction. */
  struct Ray4
  {
    /*! Default construction does nothing. */
    __forceinline Ray4() {}

    /*! Constructs a ray from origin, direction, and ray segment. Near
     *  has to be smaller than far. */
    __forceinline Ray4(const sse3f& org, const sse3f& dir, const ssef& tnear = zero, const ssef& tfar = inf, const ssef& time = zero, const ssei& mask = -1)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), id0(-1), id1(-1), mask(mask), time(time) {}

    /*! Tests if we hit something. */
    __forceinline operator sseb() const { return id0 != ssei(-1); }

  public:
    sse3f org;      //!< Ray origin
    sse3f dir;      //!< Ray direction
    ssef tnear;     //!< Start of ray segment 
    ssef tfar;      //!< End of ray segment   
    ssef u;         //!< Barycentric u coordinate of hit
    ssef v;         //!< Barycentric v coordinate of hit
    ssei id0;       //!< 1st primitive ID
    ssei id1;       //!< 2nd primitive ID
    sse3f Ng;       //!< Geometry normal
    ssei mask;      //!< used to mask out objects during traversal
    ssef time;      //!< Time of this ray for motion blur.
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Ray4& ray) {
    return cout << "{ " << 
      "org = " << ray.org << ", dir = " << ray.dir << ", near = " << ray.tnear << ", far = " << ray.tfar << ", time = " << ray.time << ", " <<
      "id0 = " << ray.id0 << ", id1 = " << ray.id1 <<  ", " << "u = " << ray.u <<  ", v = " << ray.v << ", Ng = " << ray.Ng << " }";
  }
}

#endif
