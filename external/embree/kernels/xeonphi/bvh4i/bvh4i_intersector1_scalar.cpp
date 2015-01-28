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

#include "bvh4i_intersector1_scalar.h"
#include "geometry/triangle1.h"
#include "geometry/virtual_accel_intersector1.h"

namespace embree
{
  namespace isa
  {    
    void BVH4iIntersector1Scalar::intersect(BVH4i* bvh, Ray& ray)
    {
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const Vec3fa rdir  = rcp_safe(ray.dir);
      stack_dist[0] = pos_inf;
      stack_dist[1] = pos_inf;

      const Node      * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;      
      stack_node[1] = bvh->root;

      size_t sindex = 2;

      const Vec3fa org_xyz      = ray.org;
      const Vec3fa dir_xyz      = ray.dir;
      const Vec3fa rdir_xyz     = rdir;
      const Vec3fa org_rdir_xyz = org_xyz * rdir_xyz;
	  
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;
            
	  while (1) 
	    {

	      /* test if this is a leaf node */
	      if (unlikely(curNode.isLeaf(BVH4i::leaf_mask))) break;
        
	      
	      const Node* __restrict__ const node   = curNode.node(nodes);	      
	      const Vec3fa* __restrict const plower = (Vec3fa*)node->lower;
	      const Vec3fa* __restrict const pupper = (Vec3fa*)node->upper;

	      sindex--;
	      curNode = stack_node[sindex]; // early pop of next node

	      /* intersect single ray with 4 bounding boxes */
	      float tDist[4];
	      unsigned int tNode[4];

	      size_t hiti = 0;
	      for (size_t i=0;i<4;i++)
		{
		  const Vec3fa tLowerXYZ = plower[i] * rdir_xyz - org_rdir_xyz;
		  const Vec3fa tUpperXYZ = pupper[i] * rdir_xyz - org_rdir_xyz;
		  const Vec3fa tLower = min(tLowerXYZ,tUpperXYZ);
		  const Vec3fa tUpper = max(tLowerXYZ,tUpperXYZ);
		  const float tnear = max(max(tLower.x,tLower.y),max(tLower.z,ray.tnear));
		  const float tfar  = min(min(tUpper.x,tUpper.y),min(tUpper.z,ray.tfar));
		  tDist[i] = tnear;
		  tNode[i] = node->lower[i].child;
		  if (tnear <= tfar) hiti |= (unsigned int)1 << i;
		}

	      /* if no child is hit, continue with early popped child */
	      if (unlikely(hiti == 0)) continue;
	      sindex++;
        
	      const size_t pos_first = bitscan64(hiti);
	      const size_t num_hitm = countbits(hiti); 
        
	      /* if a single child is hit, continue with that child */
	      curNode = tNode[pos_first];
	      if (likely(num_hitm == 1)) continue;
        
	      /* if two children are hit, push in correct order */
	      const size_t pos_second = bitscan64(pos_first,hiti);
	      if (likely(num_hitm == 2))
		{
		  const float dist_first  = tDist[pos_first];
		  const float dist_second = tDist[pos_second];
		  const unsigned int node_first  = curNode;
		  const unsigned int node_second = tNode[pos_second];
          
		  if (dist_first <= dist_second)
		    {
		      stack_node[sindex] = node_second;
		      stack_dist[sindex] = dist_second;                      
		      sindex++;
		      assert(sindex < 3*BVH4i::maxDepth+1);
		      continue;
		    }
		  else
		    {
		      stack_node[sindex] = curNode;
		      stack_dist[sindex] = dist_first;
		      curNode = node_second;
		      sindex++;
		      assert(sindex < 3*BVH4i::maxDepth+1);
		      continue;
		    }
		}

	      float curDist = tDist[pos_first];
	      long index = pos_first;
	      while((index = bitscan64(index,hiti)) != BITSCAN_NO_BIT_SET_64)	    
		{
		  if (curDist <= tDist[index])
		    {
		      stack_node[sindex] = tNode[index];
		      stack_dist[sindex] = tDist[index];
		    }
		  else
		    {
		      stack_node[sindex] = curNode;
		      stack_dist[sindex] = curDist;
		      
		      curNode = tNode[index];
		      curDist = tDist[index];
		    }
		  sindex++;
		  assert(sindex < 3*BVH4i::maxDepth+1);
		}        
	    }
	  
	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) 
	    {
	      break;
	    }

	  const Triangle1* tptr  = (Triangle1*) curNode.leaf(accel);

	  for (size_t i=0;i<4;i++)
	    {
	      const Vec3fa tri_v0 = tptr[i].v0;
	      const Vec3fa tri_v1 = tptr[i].v1;
	      const Vec3fa tri_v2 = tptr[i].v2;
	      const Vec3fa tri_Ng = tptr[i].Ng;

	      /* calculate denominator */
	      const Vec3fa O = ray.org;
	      const Vec3fa D = ray.dir;
	      const Vec3fa C = tri_v0 - O;
	      const Vec3fa R = cross(D,C);
	      const float den = dot(tri_Ng,D);
	      const float absDen = abs(den);
	      const float sgnDen = den < 0.0f ? -1.0f : 1.0f;
	      const Vec3fa e1 = tri_v0-tri_v1;
	      const Vec3fa e2 = tri_v2-tri_v0;

	      /* perform edge tests */
	      const float U = dot(R,e2) * sgnDen;
	      if (unlikely(U < 0.0f)) continue;
	      const float V = dot(R,e1) * sgnDen;
	      if (unlikely(V < 0.0f)) continue;
	      const float W = absDen-U-V;
	      if (unlikely(W < 0.0f)) continue;
      
	      /* perform depth test */
	      const float T = dot(tri_Ng,C) * sgnDen;
	      if (unlikely(absDen*ray.tfar < T)) continue;
	      if (unlikely(T < absDen*ray.tnear)) continue;

	      /* update hit information */
	      const float rcpAbsDen = rcp(absDen);
	      ray.u   = U * rcpAbsDen;
	      ray.v   = V * rcpAbsDen;
	      ray.tfar = T * rcpAbsDen;
	      ray.Ng  = tri_Ng;
	      ray.geomID = tptr[i].geomID();
	      ray.primID = tptr[i].primID();
	      
	    }
	}
    }


    void BVH4iIntersector1Scalar::occluded(BVH4i* bvh, Ray& ray)
    {

      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const Vec3fa rdir  = rcp_safe(ray.dir);
      const Node      * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;      
      stack_node[1] = bvh->root;

      size_t sindex = 2;

      const Vec3fa org_xyz      = ray.org;
      const Vec3fa dir_xyz      = ray.dir;
      const Vec3fa rdir_xyz     = rdir;
      const Vec3fa org_rdir_xyz = org_xyz * rdir_xyz;
	  
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;
            
	  while (1) 
	    {

	      /* test if this is a leaf node */
	      if (unlikely(curNode.isLeaf(BVH4i::leaf_mask))) break;
        
	      
	      const Node* __restrict__ const node   = curNode.node(nodes);	      
	      const Vec3fa* __restrict const plower = (Vec3fa*)node->lower;
	      const Vec3fa* __restrict const pupper = (Vec3fa*)node->upper;

	      sindex--;
	      curNode = stack_node[sindex]; // early pop of next node

	      /* intersect single ray with 4 bounding boxes */
	      float tDist[4];
	      unsigned int tNode[4];

	      size_t hiti = 0;
	      for (size_t i=0;i<4;i++)
		{
		  const Vec3fa tLowerXYZ = plower[i] * rdir_xyz - org_rdir_xyz;
		  const Vec3fa tUpperXYZ = pupper[i] * rdir_xyz - org_rdir_xyz;
		  const Vec3fa tLower = min(tLowerXYZ,tUpperXYZ);
		  const Vec3fa tUpper = max(tLowerXYZ,tUpperXYZ);
		  const float tnear = max(max(tLower.x,tLower.y),max(tLower.z,ray.tnear));
		  const float tfar  = min(min(tUpper.x,tUpper.y),min(tUpper.z,ray.tfar));
		  tDist[i] = tnear;
		  tNode[i] = node->lower[i].child;
		  if (tnear <= tfar) hiti |= (unsigned int)1 << i;
		}

	      /* if no child is hit, continue with early popped child */
	      if (unlikely(hiti == 0)) continue;
	      sindex++;
        
	      const size_t pos_first = bitscan64(hiti);
	      const size_t num_hitm = countbits(hiti); 
        
	      /* if a single child is hit, continue with that child */
	      curNode = tNode[pos_first];
	      if (likely(num_hitm == 1)) continue;
        
	      /* if two children are hit, push in correct order */
	      const size_t pos_second = bitscan64(pos_first,hiti);
	      if (likely(num_hitm == 2))
		{
		  const float dist_first  = tDist[pos_first];
		  const float dist_second = tDist[pos_second];
		  const unsigned int node_first  = curNode;
		  const unsigned int node_second = tNode[pos_second];
          
		  if (dist_first <= dist_second)
		    {
		      stack_node[sindex] = node_second;
		      sindex++;
		      assert(sindex < 3*BVH4i::maxDepth+1);
		      continue;
		    }
		  else
		    {
		      stack_node[sindex] = curNode;
		      curNode = node_second;
		      sindex++;
		      assert(sindex < 3*BVH4i::maxDepth+1);
		      continue;
		    }
		}

	      float curDist = tDist[pos_first];
	      long index = pos_first;
	      while((index = bitscan64(index,hiti)) != BITSCAN_NO_BIT_SET_64)	    
		{
		  if (curDist <= tDist[index])
		    {
		      stack_node[sindex] = tNode[index];
		    }
		  else
		    {
		      stack_node[sindex] = curNode;		      
		      curNode = tNode[index];
		      curDist = tDist[index];
		    }
		  sindex++;
		  assert(sindex < 3*BVH4i::maxDepth+1);
		}        
	    }
	  
	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) 
	    {
	      break;
	    }

	  const Triangle1* tptr  = (Triangle1*) curNode.leaf(accel);

	  for (size_t i=0;i<4;i++)
	    {
	      const Vec3fa tri_v0 = tptr[i].v0;
	      const Vec3fa tri_v1 = tptr[i].v1;
	      const Vec3fa tri_v2 = tptr[i].v2;
	      const Vec3fa tri_Ng = tptr[i].Ng;

	      /* calculate denominator */
	      const Vec3fa O = ray.org;
	      const Vec3fa D = ray.dir;
	      const Vec3fa C = tri_v0 - O;
	      const Vec3fa R = cross(D,C);
	      const float den = dot(tri_Ng,D);
	      const float absDen = abs(den);
	      const float sgnDen = den < 0.0f ? -1.0f : 1.0f;
	      const Vec3fa e1 = tri_v0-tri_v1;
	      const Vec3fa e2 = tri_v2-tri_v0;

	      /* perform edge tests */
	      const float U = dot(R,e2) * sgnDen;
	      if (unlikely(U < 0.0f)) continue;
	      const float V = dot(R,e1) * sgnDen;
	      if (unlikely(V < 0.0f)) continue;
	      const float W = absDen-U-V;
	      if (unlikely(W < 0.0f)) continue;
      
	      /* perform depth test */
	      const float T = dot(tri_Ng,C) * sgnDen;
	      if (unlikely(absDen*ray.tfar < T)) continue;
	      if (unlikely(T < absDen*ray.tnear)) continue;
	      
	      ray.geomID = 0;
	      return;
	    }
	}
    }


    DEFINE_INTERSECTOR1    (BVH4iTriangle1Intersector1Scalar, BVH4iIntersector1Scalar);
    DEFINE_INTERSECTOR1    (BVH4iVirtualIntersector1Scalar, BVH4iIntersector1Scalar);

  }
}
