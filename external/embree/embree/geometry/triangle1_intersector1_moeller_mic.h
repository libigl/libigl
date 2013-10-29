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

#ifndef __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR1_KNC2_H__
#define __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR1_KNC2_H__

#include "triangle1.h"
#include "../common/ray.h"
#include "../common/ray16.h"

namespace embree
{
  /*! Moeller Trumbore ray/triangle intersector for a single ray with
   *  individual pre-calculated triangles. */
  struct Triangle1Intersector1MoellerTrumboreMIC
  {
    typedef Triangle1 Triangle;

    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersects a ray with 4 triangles and updates the hit. */
    static __forceinline bool intersect(Ray& ray, const Triangle1* triangle, size_t items, const Vec3fa* vertices)
    {
      /* prefetch triangle */
      assert(items > 0);
      assert(items <= 4);
      prefetch<PFHINT_L1>(triangle + 3);
      prefetch<PFHINT_L1>(triangle + 2);
      prefetch<PFHINT_L1>(triangle + 1);
      prefetch<PFHINT_L1>(triangle + 0);

      /* load ray */
      const mic_f org_xyz      = loadAOS(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z);
      mic_f       min_dist_xyz = upconv1f(&ray.tnear);
      mic_f       max_dist_xyz = upconv1f(&ray.tfar);

      /* calculate edges, geometry normal, and determinant */
      const mic_i and_mask = mic_i::zlc4();
      const mic_f p0 = gather_4f_zlc(and_mask,(float*)&triangle[0].v0,(float*)&triangle[1].v0,(float*)&triangle[2].v0,(float*)&triangle[3].v0);
      const mic_f e1 = gather_4f_zlc(and_mask,(float*)&triangle[0].e1,(float*)&triangle[1].e1,(float*)&triangle[2].e1,(float*)&triangle[3].e1);
      const mic_f e2 = gather_4f_zlc(and_mask,(float*)&triangle[0].e2,(float*)&triangle[1].e2,(float*)&triangle[2].e2,(float*)&triangle[3].e2);
      const mic_f Ng = gather_4f_zlc(and_mask,(float*)&triangle[0].Ng,(float*)&triangle[1].Ng,(float*)&triangle[2].Ng,(float*)&triangle[3].Ng);
      //const mic_f Ng = lcross_xyz(e1,e2);
      
      /* calculate determinant */
      const mic_f C = p0 - org_xyz;
      const mic_f R = lcross_xyz(dir_xyz,C);
      const mic_f det = ldot3_xyz(dir_xyz,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(0x1111,det);
            
      /* perform edge tests */
      const mic_f u = ldot3_xyz(e2,R)*rcpDet;
      const mic_f v = ldot3_xyz(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) return false;

      /* perform depth test */
      const mic_f t = ldot3_xyz(C,Ng)*rcpDet;
      valid = lt(valid,min_dist_xyz,t);
      valid = gt(valid,max_dist_xyz,t);
      if (likely(none(valid))) return false;

      /* update hit information */
      max_dist_xyz = sel(valid,t,max_dist_xyz);
      const mic_f min_dist = set_min16(max_dist_xyz);
      valid = eq(valid,min_dist,max_dist_xyz);
      max_dist_xyz = min_dist;

      const unsigned long i = bsf64(valid);
      valid = andn(valid,valid-1);
      assert(countbits(valid) == 1);

      const mic_f gnormalx = mic_f(Ng[i+0]);
      const mic_f gnormaly = mic_f(Ng[i+1]);
      const mic_f gnormalz = mic_f(Ng[i+2]);
      
      compactustore16f_low(valid,&ray.u,u);
      compactustore16f_low(valid,&ray.v,v);
      compactustore16f_low(valid,&ray.tfar,t);

      compactustore16f_low(valid,&ray.Ng.x,gnormalx);
      compactustore16f_low(valid,&ray.Ng.y,gnormaly);
      compactustore16f_low(valid,&ray.Ng.z,gnormalz);
      
      const Triangle1* tri_ptr = (Triangle1*)((char*)triangle+16*i);
      ray.id0 = tri_ptr->e1.a;
      ray.id1 = tri_ptr->e2.a;
      return true;
    }

    /*! Intersects a ray with 4 triangles and updates the hit. */
    static __forceinline bool intersect(long idx, const mic_f& org_xyz, const mic_f& dir_xyz, const mic_f& min_dist_xyz, mic_f& max_dist_xyz, 
                                        Ray16& ray, const Triangle1* __restrict__ triangle, size_t items, const Vec3fa* __restrict__ vertices)
    {
      prefetch<PFHINT_L1>(triangle + 3);
      prefetch<PFHINT_L1>(triangle + 2);
      prefetch<PFHINT_L1>(triangle + 1);
      prefetch<PFHINT_L1>(triangle + 0);
      assert(items > 0);
      assert(items <= 4);
      
      /* calculate edges, geometry normal, and determinant */
      const mic_i and_mask = mic_i::zlc4();
      const mic_f p0 = gather_4f_zlc(and_mask,(float*)&triangle[0].v0,(float*)&triangle[1].v0,(float*)&triangle[2].v0,(float*)&triangle[3].v0);
      const mic_f e1 = gather_4f_zlc(and_mask,(float*)&triangle[0].e1,(float*)&triangle[1].e1,(float*)&triangle[2].e1,(float*)&triangle[3].e1);
      const mic_f e2 = gather_4f_zlc(and_mask,(float*)&triangle[0].e2,(float*)&triangle[1].e2,(float*)&triangle[2].e2,(float*)&triangle[3].e2);
      const mic_f Ng = gather_4f_zlc(and_mask,(float*)&triangle[0].Ng,(float*)&triangle[1].Ng,(float*)&triangle[2].Ng,(float*)&triangle[3].Ng);
      //const mic_f Ng = lcross_xyz(e1,e2);

      /* calculate determinant */
      const mic_f C = p0 - org_xyz;
      const mic_f R = lcross_xyz(dir_xyz,C);
      const mic_f det = ldot3_xyz(dir_xyz,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(0x1111,det);
            
      /* perform edge tests */
      const mic_f u = ldot3_xyz(e2,R)*rcpDet;
      const mic_f v = ldot3_xyz(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) return false;

      /* perform depth test */
      const mic_f t = ldot3_xyz(C,Ng)*rcpDet;
      valid = lt(valid,min_dist_xyz,t);
      valid = gt(valid,max_dist_xyz,t);
      if (likely(none(valid))) return false;

      prefetch<PFHINT_L1EX>(&ray.tfar);
      prefetch<PFHINT_L1EX>(&ray.id0);
      prefetch<PFHINT_L1EX>(&ray.id1);
      prefetch<PFHINT_L1EX>(&ray.u);
      prefetch<PFHINT_L1EX>(&ray.v);
      prefetch<PFHINT_L1EX>(&ray.Ng.x);
      prefetch<PFHINT_L1EX>(&ray.Ng.y);
      prefetch<PFHINT_L1EX>(&ray.Ng.z);

      /* update hit information */
      max_dist_xyz = sel(valid,t,max_dist_xyz);
      const mic_f min_dist = set_min16(max_dist_xyz);
      valid = eq(valid,min_dist,max_dist_xyz);
      
      const unsigned long i = bsf64(valid);
      valid = andn(valid,valid-1);
      assert(countbits(valid) == 1);

      const mic_f Ngx = Ng[i+0];
      const mic_f Ngy = Ng[i+1];
      const mic_f Ngz = Ng[i+2];
      
      compactustore16f_low(valid,&ray.u[idx],u);
      compactustore16f_low(valid,&ray.v[idx],v);
      compactustore16f_low(valid,&ray.tfar[idx],t);
      compactustore16f_low(valid,&ray.Ng.x[idx],Ngx);
      compactustore16f_low(valid,&ray.Ng.y[idx],Ngy);
      compactustore16f_low(valid,&ray.Ng.z[idx],Ngz);
      
      const Triangle1* tri_ptr = (Triangle1*)((char*)triangle+(i<<4));
      ray.id0[idx] = tri_ptr->e1.a;
      ray.id1[idx] = tri_ptr->e2.a;
      return true;
    }

    /*! Tests if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle1* __restrict__ triangle, size_t items, const Vec3fa* __restrict__ vertices)
    {
      /* prefetch triangle */
      assert(items > 0);
      assert(items <= 4);
      prefetch<PFHINT_L1>(triangle + 3);
      prefetch<PFHINT_L1>(triangle + 2);
      prefetch<PFHINT_L1>(triangle + 1);
      prefetch<PFHINT_L1>(triangle + 0);
      
      /* load ray */
      const mic_f org_xyz      = loadAOS(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f min_dist_xyz = upconv1f(&ray.tnear);
      const mic_f max_dist_xyz = upconv1f(&ray.tfar);
      
      /* calculate edges, geometry normal, and determinant */
      const mic_i and_mask = mic_i::zlc4();
      const mic_f p0 = gather_4f_zlc(and_mask,(float*)&triangle[0].v0,(float*)&triangle[1].v0,(float*)&triangle[2].v0,(float*)&triangle[3].v0);
      const mic_f e1 = gather_4f_zlc(and_mask,(float*)&triangle[0].e1,(float*)&triangle[1].e1,(float*)&triangle[2].e1,(float*)&triangle[3].e1);
      const mic_f e2 = gather_4f_zlc(and_mask,(float*)&triangle[0].e2,(float*)&triangle[1].e2,(float*)&triangle[2].e2,(float*)&triangle[3].e2);
      const mic_f Ng = gather_4f_zlc(and_mask,(float*)&triangle[0].Ng,(float*)&triangle[1].Ng,(float*)&triangle[2].Ng,(float*)&triangle[3].Ng);
      //const mic_f Ng = lcross_xyz(e1,e2);
      
      /* calculate determinant */
      const mic_f C = p0 - org_xyz;
      const mic_f R = lcross_xyz(dir_xyz,C);
      const mic_f det = ldot3_xyz(dir_xyz,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(0x1111,det);
            
      /* perform edge tests */
      const mic_f u = ldot3_xyz(e2,R)*rcpDet;
      const mic_f v = ldot3_xyz(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) return false;

      /* perform depth test */
      const mic_f t = ldot3_xyz(C,Ng)*rcpDet;
      valid = lt(valid,min_dist_xyz,t);
      valid = gt(valid,max_dist_xyz,t);
      return any(valid);
    }
  };
}

#endif


