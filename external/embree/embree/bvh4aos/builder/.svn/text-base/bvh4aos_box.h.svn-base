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

#ifndef __EMBREE_BVH4AOS_MIC_BOX_H__
#define __EMBREE_BVH4AOS_MIC_BOX_H__

#include "bvh4aos_globals.h"

/* ---------------------------------- */
/* --- Extra data for QBVH nodes  --- */
/* ---------------------------------- */

namespace embree
{

  class AABBExtData 
  {
  public:
    union {
      class { 
      public:
	unsigned short items;
	unsigned short flags;
      };
      union {
	unsigned int children;
	unsigned int t;
	int v;
      };
    };
    AABBExtData() {}
  };

  /* ------------ */
  /* --- AABB --- */
  /* ------------ */

  __ALIGN(32)
    class AABB
    {
    public:
      float m_min[3];
      AABBExtData ext_min;
      float m_max[3];
      AABBExtData ext_max;

      _INLINE float area() const
      {
	const float x = m_max[0] - m_min[0];
	const float y = m_max[1] - m_min[1];
	const float z = m_max[2] - m_min[2];
	return 2.0f * (x*y+x*z+y*z); 
      }

      _INLINE float centroid(const int dim) const
      {
	return (m_min[dim] + m_max[dim]) * 0.5f;
      }

      _INLINE float centroid_2(const int dim) const
      {
	return m_min[dim] + m_max[dim];
      }

      _INLINE void set(const Vec3fa &v)
      {
	m_min[0] = v.x;
	m_max[0] = v.x;
	m_min[1] = v.y;
	m_max[1] = v.y;
	m_min[2] = v.z;
	m_max[2] = v.z;    
      }

      _INLINE void setEmpty()
      {
	m_min[0] = float_infinity;
	m_min[1] = float_infinity;
	m_min[2] = float_infinity;
	ext_min.t = 0;
	m_max[0] = float_minus_infinity;
	m_max[1] = float_minus_infinity;
	m_max[2] = float_minus_infinity;
	ext_max.t = 0;
      }

      _INLINE void resetBounds()
      {
	m_min[0] = float_infinity;
	m_min[1] = float_infinity;
	m_min[2] = float_infinity;
	m_max[0] = float_minus_infinity;
	m_max[1] = float_minus_infinity;
	m_max[2] = float_minus_infinity;
      }


      _INLINE void setZero()
      {
	m_min[0] = 0.0f;
	m_min[1] = 0.0f; 
	m_min[2] = 0.0f; 
	ext_min.t = 0;
	m_max[0] = 0.0f; 
	m_max[1] = 0.0f; 
	m_max[2] = 0.0f; 
	ext_max.t = 0;
      }

      _INLINE void setToInfPoint()
      {
	m_min[0] = 1E14;
	m_min[1] = 1E14;
	m_min[2] = 1E14;
	ext_min.t = 0;
	m_max[0] = 1E14;
	m_max[1] = 1E14;
	m_max[2] = 1E14;
	ext_max.t = 0;
      }

      _INLINE void extend(const Vec3fa &v)
      {
	m_min[0] = min(m_min[0],v.x);
	m_max[0] = max(m_max[0],v.x);
	m_min[1] = min(m_min[1],v.y);
	m_max[1] = max(m_max[1],v.y);
	m_min[2] = min(m_min[2],v.z);
	m_max[2] = max(m_max[2],v.z);
      }

      _INLINE void extend(const AABB &aabb)
      {
	m_min[0] = min(m_min[0],aabb.m_min[0]);
	m_max[0] = max(m_max[0],aabb.m_max[0]);
	m_min[1] = min(m_min[1],aabb.m_min[1]);
	m_max[1] = max(m_max[1],aabb.m_max[1]);
	m_min[2] = min(m_min[2],aabb.m_min[2]);
	m_max[2] = max(m_max[2],aabb.m_max[2]);
      }

      _INLINE bool enclose(const AABB &aabb) const
      {
	return
	  m_min[0] <= aabb.m_min[0] && aabb.m_max[0] <= m_max[0] &&
	  m_min[1] <= aabb.m_min[1] && aabb.m_max[1] <= m_max[1] &&
	  m_min[2] <= aabb.m_min[2] && aabb.m_max[2] <= m_max[2];      
      }

      _INLINE bool enclose(const Vec3fa &v) const
      {
	return
	  m_min[0] <= v[0] && v[0] <= m_max[0] &&
	  m_min[1] <= v[1] && v[1] <= m_max[1] &&
	  m_min[2] <= v[2] && v[2] <= m_max[2];      
      }

      _INLINE bool isEmpty()
      {
	if (m_min[0] == float_infinity &&
	    m_min[1] == float_infinity &&
	    m_min[2] == float_infinity &&
	    m_max[0] == float_minus_infinity &&
	    m_max[1] == float_minus_infinity &&
	    m_max[2] == float_minus_infinity)
	  return true;
	return false;
      }

      _INLINE void operator=(const AABB& v) { 
	const mic_f b = uload16f_low(v.m_min);
	compactustore16f_low(0xff,m_min,b);
      }

    };

  _INLINE void xchg(AABB& a, AABB& b)
  {
    const mic_f aa = uload16f_low(a.m_min);
    const mic_f bb = uload16f_low(b.m_min);
    compactustore16f_low(0xff,b.m_min,aa);
    compactustore16f_low(0xff,a.m_min,bb);
  }

  _INLINE std::ostream &operator<<(std::ostream &o, const AABB &v)
  {
    o << "min [" << v.m_min[0] << ", " << v.m_min[1] << ", " << v.m_min[2] << ", " << v.ext_min.t <<"] ";
    o << "max [" << v.m_max[0] << ", " << v.m_max[1] << ", " << v.m_max[2] << ", " << v.ext_max.t <<"] ";
    return o;
  } 

  _INLINE bool equal(AABB &a, AABB &b)
  {
    if (a.m_min[0] == b.m_min[0] &&
	a.m_min[1] == b.m_min[1] &&
	a.m_min[2] == b.m_min[2] &&
	a.m_max[0] == b.m_max[0] &&
	a.m_max[1] == b.m_max[1] &&
	a.m_max[2] == b.m_max[2])
      return true;
    return false;
  };

};

#endif
