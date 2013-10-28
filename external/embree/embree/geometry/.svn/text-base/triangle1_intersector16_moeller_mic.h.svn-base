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

#ifndef __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR16_H__
#define __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR16_H__

#include "triangle1.h"
#include "triangle1_intersector1_moeller_mic.h"

#include "../common/ray16.h"

namespace embree
{
  /*! Moeller Trumbore ray/triangle intersector for a single ray with
   *  individual pre-gathered triangles. */
  struct Triangle1Intersector16MoellerTrumboreMIC
  {
    typedef Triangle1 Triangle;
    typedef Triangle1Intersector1MoellerTrumboreMIC TriangleIntersector1; 

    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersects 16 rays with 4 triangles and updates the hit. */
    static __forceinline void intersect(const mic_m valid_i, Ray16& ray, const Triangle1* triangle, const Vec3fa* vertices)
    {
      /* load ray into registers */
      const mic3f org(ray.org.x,ray.org.y,ray.org.z);
      const mic3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f min_distance = ray.tnear;
      const mic_f max_distance = ray.tfar;

      /* calculate edges and geometry normal */
      const mic_f _p0 = upconv4f((float*)&triangle->v0);
      const mic_f _e1 = upconv4f((float*)&triangle->e1);
      const mic_f _e2 = upconv4f((float*)&triangle->e2);
      const mic_f _Ng = upconv4f((float*)&triangle->Ng);
      //const mic_f _Ng = lcross_xyz(_e1,_e2);
      
      /* change data layout */
      const mic3f p0(swAAAA(_p0),swBBBB(_p0),swCCCC(_p0));
      const mic3f e1(swAAAA(_e1),swBBBB(_e1),swCCCC(_e1));
      const mic3f e2(swAAAA(_e2),swBBBB(_e2),swCCCC(_e2));
      const mic3f Ng(swAAAA(_Ng),swBBBB(_Ng),swCCCC(_Ng));
      
      /* calculate determinant */
      const mic3f C = p0 - org;
      const mic3f R = cross(dir,C);
      const mic_f det = dot(dir,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(valid_i,det);
      
      /* perform edge tests */
      const mic_f u = dot(e2,R)*rcpDet;
      const mic_f v = dot(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const mic_f t = dot(C,Ng)*rcpDet;
      valid = lt(valid,min_distance,t);
      valid = gt(valid,max_distance,t);
      const mic_i id0 = upconv1i((const int *)&triangle->e1.a);
      const mic_i id1 = upconv1i((const int *)&triangle->e2.a);
      if (likely(none(valid))) return;
      
      /* update hit information */
      prefetch<PFHINT_L1EX>(&ray.u);
      prefetch<PFHINT_L1EX>(&ray.v);
      prefetch<PFHINT_L1EX>(&ray.tfar);
      prefetch<PFHINT_L1EX>(&ray.id0);
      prefetch<PFHINT_L1EX>(&ray.id1);
      store16f(valid,&ray.Ng.x,Ng.x);
      store16f(valid,&ray.Ng.y,Ng.y);
      store16f(valid,&ray.Ng.z,Ng.z);
      store16f(valid,&ray.u,u);
      store16f(valid,&ray.v,v);
      store16f(valid,&ray.tfar,t);
      store16i(valid,&ray.id0,id0);
      store16i(valid,&ray.id1,id1);
    }

    /*! Intersects 16 rays with 4 triangles and updates the hit. */
    static __forceinline void intersect(const mic_m valid, Ray16& ray, const Triangle1* triangle, size_t items, const Vec3fa* vertices)
    {
      assert(items > 0);
      assert(items <= 4);

      prefetch<PFHINT_L1>((mic_f*)triangle +  0); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  1); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  2); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  3); 

      for (unsigned int i=0; i<items; i++, triangle++) {
        prefetch<PFHINT_L1>((mic_f*)triangle +  1); 
        intersect(valid,ray,triangle,vertices);
      }
    }

    /*! Tests if the ray is occluded by one of the triangles. */
    static __forceinline mic_m occluded(const mic_m valid_i, const Ray16& ray, const Triangle1* triangle, const Vec3fa* vertices)
    {
      /* load ray into registers */
      const mic3f org(ray.org.x,ray.org.y,ray.org.z);
      const mic3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f min_distance = ray.tnear;
      const mic_f max_distance = ray.tfar;

      /* calculate edges and geometry normal */
      const mic_f _p0 = upconv4f((float*)&triangle->v0);
      const mic_f _e1 = upconv4f((float*)&triangle->e1);
      const mic_f _e2 = upconv4f((float*)&triangle->e2);
      const mic_f _Ng = upconv4f((float*)&triangle->Ng);
      //const mic_f _Ng = lcross_xyz(_e1,_e2);
      
      /* change data layout */
      const mic3f p0(swAAAA(_p0),swBBBB(_p0),swCCCC(_p0));
      const mic3f e1(swAAAA(_e1),swBBBB(_e1),swCCCC(_e1));
      const mic3f e2(swAAAA(_e2),swBBBB(_e2),swCCCC(_e2));
      const mic3f Ng(swAAAA(_Ng),swBBBB(_Ng),swCCCC(_Ng));
      
      /* calculate determinant */
      const mic3f C = p0 - org;
      const mic3f R = cross(dir,C);
      const mic_f det = dot(dir,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(valid_i,det);
      
      /* perform edge tests */
      const mic_f u = dot(e2,R)*rcpDet;
      const mic_f v = dot(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) return 0;
      
      /* perform depth test */
      const mic_f t = dot(C,Ng)*rcpDet;
      valid = lt(valid,min_distance,t);
      valid = gt(valid,max_distance,t);
      return valid;
    }

    /*! Tests if the ray is occluded by one of the triangles. */
    static __forceinline mic_m occluded(const mic_m valid, const Ray16& ray, const Triangle1* triangle, size_t items, const Vec3fa* vertices)
    {
      assert(items > 0);
      assert(items <= 4);
      prefetch<PFHINT_L1>((mic_f*)triangle +  0); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  1); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  2); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  3); 
      
      mic_m valid_o = 0;
      for (unsigned int i=0;i<items;i++,triangle++) 
      {
        prefetch<PFHINT_L1>((mic_f*)triangle +  1); 
        valid_o |= occluded(valid,ray,triangle,vertices);
        if (unlikely((valid & ~valid_o) == 0)) break;
      }
      return valid_o;
    }
  };
}

#endif


