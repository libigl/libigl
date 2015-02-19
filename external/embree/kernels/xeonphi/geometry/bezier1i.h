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

#include "common/default.h"
#include "primitive.h"

#define EVAL_BOUNDS 1

namespace embree
{
  extern mic4f coeff0;
  extern mic4f coeff1;
  extern mic4f coeff01;

  extern mic4f coeff_P0;
  extern mic4f coeff_P1;
  extern mic4f coeff_P2;
  extern mic4f coeff_P3;

  __forceinline mic3f convert(const LinearSpace3fa &mat)
  {
    return mic3f(broadcast4to16f((float*)&mat.vx),
		 broadcast4to16f((float*)&mat.vy),
		 broadcast4to16f((float*)&mat.vz));
  }

  __forceinline mic_f xfmPoint4f(const Vec3fa &p, 
				 const mic_f &c0,
				 const mic_f &c1,
				 const mic_f &c2)
  {
    const mic_f xyz  = c0 * mic_f(p.x) + c1 * mic_f(p.y) + c2 * mic_f(p.z);
    const mic_f xyzw = select(0x7777,xyz,mic_f(p.w));
    return xyzw;
  }


    static __forceinline mic4f eval16(const mic_f &p0,
				      const mic_f &p1,
				      const mic_f &p2,
				      const mic_f &p3)
    {
      const mic_f c0 = coeff01.x;
      const mic_f c1 = coeff01.y;
      const mic_f c2 = coeff01.z;
      const mic_f c3 = coeff01.w;

      const mic_f x = c0 * swAAAA(p0) + c1 * swAAAA(p1) + c2 * swAAAA(p2) + c3 * swAAAA(p3);
      const mic_f y = c0 * swBBBB(p0) + c1 * swBBBB(p1) + c2 * swBBBB(p2) + c3 * swBBBB(p3);
      const mic_f z = c0 * swCCCC(p0) + c1 * swCCCC(p1) + c2 * swCCCC(p2) + c3 * swCCCC(p3);
      const mic_f w = c0 * swDDDD(p0) + c1 * swDDDD(p1) + c2 * swDDDD(p2) + c3 * swDDDD(p3);
      return mic4f(x,y,z,w);
    }


  struct __aligned(16) Bezier1i
  {
  public:

    /*! Default constructor. */
    __forceinline Bezier1i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Bezier1i (const Vec3fa* p, const unsigned int geomID, const unsigned int primID)
      : p(p), geomID(geomID), primID(primID) {}

    /*! calculate the bounds of the triangle */
    __forceinline BBox3fa bounds() const {
      const BBox3fa b = merge(BBox3fa(p[0]),BBox3fa(p[1]),BBox3fa(p[2]),BBox3fa(p[3]));
      return enlarge(b,Vec3fa(b.upper.w));
    }

    __forceinline mic2f getBounds() const 
    {
      const mic_f v0 = broadcast4to16f((float*)&p[0]);
      const mic_f v1 = broadcast4to16f((float*)&p[1]);
      const mic_f v2 = broadcast4to16f((float*)&p[2]);
      const mic_f v3 = broadcast4to16f((float*)&p[3]);
      
#if EVAL_BOUNDS==1
      const mic4f v = eval16(v0,v1,v2,v3);
      const mic_f min_x = min(vreduce_min(v.x),v3[0]);
      const mic_f max_x = max(vreduce_max(v.x),v3[0]);
      const mic_f min_y = min(vreduce_min(v.y),v3[1]);
      const mic_f max_y = max(vreduce_max(v.y),v3[1]);
      const mic_f min_z = min(vreduce_min(v.z),v3[2]);
      const mic_f max_z = max(vreduce_max(v.z),v3[2]);
      const mic_f b_min = select(0x4444,min_z,select(0x2222,min_y,min_x));
      const mic_f b_max = select(0x4444,max_z,select(0x2222,max_y,max_x));

      const mic_f r_max = max(max(v0,v1),max(v2,v3));
      const mic_f b_min_r = b_min - swDDDD(r_max);
      const mic_f b_max_r = b_max + swDDDD(r_max);

#else
      const mic_f b_min = min(min(v0,v1),min(v2,v3));
      const mic_f b_max = max(max(v0,v1),max(v2,v3));
      
      const mic_f b_min_r = b_min - swDDDD(b_max);
      const mic_f b_max_r = b_max + swDDDD(b_max);
#endif      
      return mic2f(b_min_r,b_max_r);
    }


    __forceinline mic2f getBounds(LinearSpace3fa &xfm) const 
    {
      const Vec3fa q0 = xfmPoint(xfm,p[0]);
      const Vec3fa q1 = xfmPoint(xfm,p[1]);
      const Vec3fa q2 = xfmPoint(xfm,p[2]);
      const Vec3fa q3 = xfmPoint(xfm,p[3]);
      
      const Vec3fa b_min = min(min(q0,q1),min(q2,q3));
      const Vec3fa b_max = max(max(q0,q1),max(q2,q3));

      const Vec3fa max_radius = max(max(p[0].w,p[1].w),max(p[2].w,p[3].w));
      
      const Vec3fa b_min_r = b_min - max_radius;
      const Vec3fa b_max_r = b_max + max_radius;

      const mic_f b_lower = broadcast4to16f((float*)&b_min_r);
      const mic_f b_upper = broadcast4to16f((float*)&b_max_r);
      
      return mic2f(b_lower,b_upper);
    }

    __forceinline mic2f getBounds(const mic_f &c0,const mic_f &c1,const mic_f &c2) const 
    {
      const mic_f v0 = xfmPoint4f(p[0],c0,c1,c2);
      const mic_f v1 = xfmPoint4f(p[1],c0,c1,c2);
      const mic_f v2 = xfmPoint4f(p[2],c0,c1,c2);
      const mic_f v3 = xfmPoint4f(p[3],c0,c1,c2);

#if EVAL_BOUNDS==1
      const mic4f v = eval16(v0,v1,v2,v3);
      const mic_f min_x = min(vreduce_min(v.x),v3[0]);
      const mic_f max_x = max(vreduce_max(v.x),v3[0]);
      const mic_f min_y = min(vreduce_min(v.y),v3[1]);
      const mic_f max_y = max(vreduce_max(v.y),v3[1]);
      const mic_f min_z = min(vreduce_min(v.z),v3[2]);
      const mic_f max_z = max(vreduce_max(v.z),v3[2]);
      const mic_f b_min = select(0x4444,min_z,select(0x2222,min_y,min_x));
      const mic_f b_max = select(0x4444,max_z,select(0x2222,max_y,max_x));

      const mic_f r_max = max(max(v0,v1),max(v2,v3));
      const mic_f b_min_r = b_min - swDDDD(r_max);
      const mic_f b_max_r = b_max + swDDDD(r_max);

#else
      const mic_f b_min = min(min(v0,v1),min(v2,v3));
      const mic_f b_max = max(max(v0,v1),max(v2,v3));
      const mic_f b_min_r = b_min - swDDDD(b_max);
      const mic_f b_max_r = b_max + swDDDD(b_max);
#endif      
      
      return mic2f(b_min_r,b_max_r);
    }       

    template<int HINT>
    __forceinline void prefetchControlPoints() const {
      prefetch<HINT>(p + 0);
      prefetch<HINT>(p + 3);
    }

    const Vec3fa* p;      //!< pointer to first control point (x,y,z,r)
    unsigned int geomID;  //!< geometry ID
    unsigned int primID;  //!< primitive ID
  };

};
