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
      : org(org), dir(dir), tnear(tnear), tfar(tfar), geomID(-1), primID(-1), instID(-1), mask(mask), time(time)  {}

    /*! Tests if we hit something. */
    __forceinline operator avxb() const { return geomID != avxi(-1); }

  public:
    avx3f org;      //!< Ray origin
    avx3f dir;      //!< Ray direction
    avxf tnear;     //!< Start of ray segment 
    avxf tfar;      //!< End of ray segment   
    avxf time;      //!< Time of this ray for motion blur.
    avxi mask;      //!< used to mask out objects during traversal

  public:
    avx3f Ng;       //!< Geometry normal
    avxf u;         //!< Barycentric u coordinate of hit
    avxf v;         //!< Barycentric v coordinate of hit
    avxi geomID;    //!< geometry ID
    avxi primID;    //!< primitive ID
    avxi instID;    //!< instance ID
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Ray8& ray) {
    return cout << "{ " << 
      "org = " << ray.org << ", dir = " << ray.dir << ", near = " << ray.tnear << ", far = " << ray.tfar << ", time = " << ray.time << ", " <<
      "instID = " << ray.instID << ", geomID = " << ray.geomID << ", primID = " << ray.primID <<  ", " << "u = " << ray.u <<  ", v = " << ray.v << ", Ng = " << ray.Ng << " }";
  }
}

#endif
