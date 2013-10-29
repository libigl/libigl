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

#ifndef __EMBREE_RAY8_H__
#define __EMBREE_RAY8_H__

#include "ray.h"

namespace embree
{
  /*! Ray structure. Contains all information of 8 rays including
   *  precomputed reciprocal direction. */
  struct Ray8
  {
    /*! Default construction does nothing. */
    __forceinline Ray8() {}

    /*! Constructs a ray from origin, direction, and ray segment. Near
     *  has to be smaller than far. */
    __forceinline Ray8(const avx3f& org, const avx3f& dir, const avxf& tnear = zero, const avxf& tfar = inf, const avxf& time = zero, const avxi& mask = -1)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), id0(-1), id1(-1), mask(mask), time(time)  {}

    /*! Tests if we hit something. */
    __forceinline operator avxb() const { return id0 != avxi(-1); }

  public:
    avx3f org;      //!< Ray origin
    avx3f dir;      //!< Ray direction
    avxf tnear;     //!< Start of ray segment 
    avxf tfar;      //!< End of ray segment   
    avxf u;         //!< Barycentric u coordinate of hit
    avxf v;         //!< Barycentric v coordinate of hit
    avxi id0;       //!< 1st primitive ID
    avxi id1;       //!< 2nd primitive ID
    avx3f Ng;       //!< Geometry normal
    avxi mask;      //!< used to mask out objects during traversal
    avxf time;      //!< Time of this ray for motion blur.
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Ray8& ray) {
    return cout << "{ " << 
      "org = " << ray.org << ", dir = " << ray.dir << ", near = " << ray.tnear << ", far = " << ray.tfar << ", time = " << ray.time << ", " <<
      "id0 = " << ray.id0 << ", id1 = " << ray.id1 <<  ", " << "u = " << ray.u <<  ", v = " << ray.v << ", Ng = " << ray.Ng << " }";
  }
}

#endif
