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

#ifndef __EMBREE_RAY16_H__
#define __EMBREE_RAY16_H__

#include "ray.h"

namespace embree
{
  /*! Ray structure. Contains all information of 16 rays including
   *  precomputed reciprocal direction. */
  struct Ray16
  {
    /*! Default construction does nothing. */
    __forceinline Ray16() {}

    /*! Constructs a ray from origin, direction, and ray segment. Near
     *  has to be smaller than far. */
    __forceinline Ray16(const mic3f& org, const mic3f& dir, const mic_f& tnear = zero, const mic_f& tfar = inf, const mic_f& time = zero, const mic_i& mask = -1)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), geomID(-1), primID(-1), mask(mask), time(time) {}

    /*! Tests if we hit something. */
    __forceinline operator mic_m() const { return geomID != mic_i(-1); }

  public:
    mic3f org;      //!< Ray origin
    mic3f dir;      //!< Ray direction
    mic_f tnear;    //!< Start of ray segment 
    mic_f tfar;     //!< End of ray segment   
    mic_f time;     //!< Time of this ray for motion blur.
    mic_i mask;     //!< used to mask out objects during traversal

  public:
    mic3f Ng;       //!< Geometry normal
    mic_f u;        //!< Barycentric u coordinate of hit
    mic_f v;        //!< Barycentric v coordinate of hit
    mic_i geomID;   //!< geometry ID
    mic_i primID;   //!< primitive ID
    mic_i instID;   //!< instance ID
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Ray16& ray) {
    return cout << "{ " << std::endl << "org = " << ray.org << std::endl << " dir = " << ray.dir << std::endl << " near = " << ray.tnear << std::endl << " far = " << ray.tfar << std::endl << " time = " << ray.time << std::endl << " " <<
      "instID = " << ray.instID << std::endl << " geomID = " << ray.geomID << std::endl << " primID = " << ray.primID <<  std::endl << " " << "u = " << ray.u <<  std::endl << " v = " << ray.v << std::endl << " Ng = " << ray.Ng << " }";
  }
}

#endif
