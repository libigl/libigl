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

#include "bvh4aos_triangle1_intersector16_single.h"
#include "bvh4aos_triangle1_intersector1.h"
#include "../geometry/triangles.h"

namespace embree
{
  void BVH4AOSTriangle1Intersector16Single::intersect(const BVH4AOSTriangle1Intersector16Single* This, Ray16& ray, const __mmask valid_i)
  {
    /* pointers to node array and triangle array */
    const mic_m valid = valid_i;
    const BVH4AOS* bvh = This->bvh;
    const Node*     const nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle1 *const triangles = (const Triangle1*) bvh->triPtr();

    long rayIndex = -1;
    while((rayIndex = bsf64(rayIndex,valid)) != MIC_NO_BIT_SET_64)	    
    {	    
      const mic3f ray_rdir     = rcp_safe(ray.dir);
      const mic_f org_xyz      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = SOAtoAOS_4f(rayIndex,ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = SOAtoAOS_4f(rayIndex,ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = upconv1f(&ray.tnear[rayIndex]); 
      const mic_f max_dist_xyz = upconv1f(&ray.tfar[rayIndex]); 
      BVH4AOSTriangle1Intersector1::intersect1(bvh,nodes,triangles,bvh->root,rayIndex,org_xyz,dir_xyz,rdir_xyz,org_rdir_xyz,min_dist_xyz,max_dist_xyz,ray);
    }
  }

  __mmask BVH4AOSTriangle1Intersector16Single::occluded(const BVH4AOSTriangle1Intersector16Single* This, Ray16& ray, const __mmask valid_i)
  {
    /* pointers to node array and triangle array */
    const mic_m valid = valid_i;
    const BVH4AOS* bvh = This->bvh;
    const Node*     const nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle1 *const triangles = (const Triangle1*) bvh->triPtr();

    mic_i not_occluded = mic_i::minus_one();
    long rayIndex = -1;
    while((rayIndex = bsf64(rayIndex,valid)) != MIC_NO_BIT_SET_64)	    
    {	    
      // === TODO: precompute SOAtoAOS transformation, load with 4x broadcast
      const mic3f ray_rdir     = rcp_safe(ray.dir);
      const mic_f org_xyz      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = SOAtoAOS_4f(rayIndex,ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = SOAtoAOS_4f(rayIndex,ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = upconv1f(&ray.tnear[rayIndex]); 
      const mic_f max_dist_xyz = upconv1f(&ray.tfar[rayIndex]); 
      if (BVH4AOSTriangle1Intersector1::occluded1(bvh,nodes,triangles,bvh->root,rayIndex,org_xyz,dir_xyz,rdir_xyz,org_rdir_xyz,min_dist_xyz,max_dist_xyz))
        not_occluded[rayIndex] = 0;
    }
    return valid & eq(not_occluded,mic_i::zero());
  }

  void BVH4AOSTriangle1Intersector16SingleRegister () 
  {
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","single_inlined","moeller",true ,BVH4AOSTriangle1Intersector16Single::create);
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","single","moeller",true ,BVH4AOSTriangle1Intersector16Single::create);

  }
}
