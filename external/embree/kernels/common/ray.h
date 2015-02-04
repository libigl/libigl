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

#include "default.h"

namespace embree
{
  /*! Ray structure. Contains all information about a ray including
   *  precomputed reciprocal direction. */
  struct Ray
  {
    /*! Default construction does nothing. */
    __forceinline Ray() {}

    /*! Constructs a ray from origin, direction, and ray segment. Near
     *  has to be smaller than far. */
    __forceinline Ray(const Vec3fa& org, const Vec3fa& dir, float tnear = zero, float tfar = inf, float time = zero, int mask = -1)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), geomID(-1), primID(-1), instID(-1), mask(mask), time(time) {}

    /*! Tests if we hit something. */
    __forceinline operator bool() const { return geomID != -1; }

  public:
    Vec3fa org;        //!< Ray origin
    Vec3fa dir;        //!< Ray direction
    float tnear;       //!< Start of ray segment
    float tfar;        //!< End of ray segment
    float time;        //!< Time of this ray for motion blur.
    int mask;          //!< used to mask out objects during traversal

    Vec3fa Ng;         //!< Not normalized geometry normal
    float u;           //!< Barycentric u coordinate of hit
    float v;           //!< Barycentric v coordinate of hit
    int geomID;        //!< geometry ID
    int primID;        //!< primitive ID
    int instID;        //!< instance ID

#if defined(__MIC__)    
    __forceinline void update(const mic_m &m_mask,
			      const mic_f &new_t,
			      const mic_f &new_u,
			      const mic_f &new_v,
			      const mic_f &new_gnormalx,
			      const mic_f &new_gnormaly,
			      const mic_f &new_gnormalz,
			      const int new_geomID,
			      const int new_primID)
    {
      geomID = new_geomID;
      primID = new_primID;

      compactustore16f_low(m_mask,&tfar,new_t);
      compactustore16f_low(m_mask,&u,new_u); 
      compactustore16f_low(m_mask,&v,new_v); 
      compactustore16f_low(m_mask,&Ng.x,new_gnormalx); 
      compactustore16f_low(m_mask,&Ng.y,new_gnormaly); 
      compactustore16f_low(m_mask,&Ng.z,new_gnormalz);       
    }
#endif
  };

  /*! Outputs ray to stream. */
  inline std::ostream& operator<<(std::ostream& cout, const Ray& ray) {
    return cout << "{ " << 
      "org = " << ray.org << ", dir = " << ray.dir << ", near = " << ray.tnear << ", far = " << ray.tfar << ", time = " << ray.time << ", " <<
      "instID = " << ray.instID << ", geomID = " << ray.geomID << ", primID = " << ray.primID <<  ", " << "u = " << ray.u <<  ", v = " << ray.v << ", Ng = " << ray.Ng << " }";
  }
}
