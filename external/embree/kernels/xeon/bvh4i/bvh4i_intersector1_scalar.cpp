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
#include "../common/stack_item.h"

#include "geometry/triangle1_intersector1_moeller.h"
#include "geometry/triangle4_intersector1_moeller.h"
#include "geometry/triangle1v_intersector1_pluecker.h"
#include "geometry/triangle4v_intersector1_pluecker.h"
#include "geometry/virtual_accel_intersector1.h"

namespace embree
{
  namespace isa
  {

    static __forceinline void intersect_vec3f(Ray& ray, const Triangle1& tri, const void* geom)
    {
      /* load triangle */
      STAT3(normal.trav_prims,1,1,1);
      const Vec3f tri_v0 = tri.v0;
      const Vec3f tri_v1 = tri.v1;
      const Vec3f tri_v2 = tri.v2;
      const Vec3f tri_Ng = tri.Ng;

      /* calculate denominator */
      const Vec3f O = ray.org;
      const Vec3f D = ray.dir;
      const Vec3f C = tri_v0 - O;
      const Vec3f R = cross(D,C);
      const float den = dot(tri_Ng,D);
      const float absDen = abs(den);
      const float sgnDen = den < 0.0f ? -1.0f : 1.0f;
      const Vec3f e1 = tri_v0-tri_v1;
      const Vec3f e2 = tri_v2-tri_v0;

      /* perform edge tests */
      const float U = dot(R,e2) * sgnDen;
      if (unlikely(U < 0.0f)) return;
      const float V = dot(R,e1) * sgnDen;
      if (unlikely(V < 0.0f)) return;
      const float W = absDen-U-V;
      if (unlikely(W < 0.0f)) return;
      
      /* perform depth test */
      const float T = dot(tri_Ng,C) * sgnDen;
      if (unlikely(absDen*ray.tfar < T)) return;
      if (unlikely(T < absDen*ray.tnear)) return;

      /* update hit information */
      const float rcpAbsDen = rcp(absDen);
      ray.u   = U * rcpAbsDen;
      ray.v   = V * rcpAbsDen;
      ray.tfar = T * rcpAbsDen;
      ray.Ng  = tri_Ng;
      ray.geomID = tri.geomID();
      ray.primID = tri.primID();
    }


    static __forceinline bool occluded_vec3f(const Ray& ray, const Triangle1& tri, const void* geom)
    {
      /* load triangle */
      STAT3(shadow.trav_prims,1,1,1);
      const Vec3f tri_v0 = tri.v0;
      const Vec3f tri_v1 = tri.v1;
      const Vec3f tri_v2 = tri.v2;
      const Vec3f tri_Ng = tri.Ng;

      /* calculate denominator */
      const Vec3f O = Vec3fa(ray.org);
      const Vec3f D = Vec3fa(ray.dir);
      const Vec3f C = tri_v0 - O;
      const Vec3f R = cross(D,C);
      const float den = dot(tri_Ng,D);
      const float absDen = abs(den);
      const float sgnDen = den < 0.0f ? -1.0f : 1.0f;
      const Vec3fa e1 = tri_v0-tri_v1;
      const Vec3fa e2 = tri_v2-tri_v0;

      /* perform edge tests */
      const float U = dot(R,e2) * sgnDen;
      if (unlikely(U < 0.0f)) return false;
      const float V = dot(R,e1) * sgnDen;
      if (unlikely(V < 0.0f)) return false;
      const float W = absDen-U-V;
      if (unlikely(W < 0.0f)) return false;
      
      /* perform depth test */
      const float T = dot(tri_Ng,C) * sgnDen;
      if (unlikely(absDen*ray.tfar < T)) return false;
      if (unlikely(T < absDen*ray.tnear)) return false;

      return true;
    }


    template<typename TriangleIntersector>
    void BVH4iIntersector1Scalar<TriangleIntersector>::intersect(const BVH4i* bvh, Ray& ray)
    {
      /*! stack state */
      StackItem stack[1+3*BVH4i::maxDepth];  //!< stack of nodes 
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      stack[0].ptr  = bvh->root;
      stack[0].dist = neg_inf;
      
      
      /*! load the ray into SIMD registers */
      const Vec3f rdir = rcp_safe(ray.dir);
      const Vec3f org_rdir = ray.org*rdir;
      
      const void* nodePtr = bvh->nodePtr();
      const void* triPtr  = bvh->triPtr();
      
      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        NodeRef cur = NodeRef(stackPtr->ptr);
        
        /*! if popped node is too far, pop next one */
        if (unlikely(stackPtr->dist > ray.tfar))
          continue;
        
        /* downtraversal loop */
        while (true)
        {
          /*! stop if we found a leaf */
          if (unlikely(cur.isLeaf())) break;
          STAT3(normal.trav_nodes,1,1,1);
          
          /*! single ray intersection with 4 boxes */
          const Node* node = cur.node(nodePtr);

	  size_t pushed = 0;
	  for (size_t i=0;i<4;i++)
	    {
	      const float nearX = node->lower_x[i] * rdir.x - org_rdir.x;
	      const float farX  = node->upper_x[i] * rdir.x - org_rdir.x;
	      const float nearY = node->lower_y[i] * rdir.y - org_rdir.y;
	      const float farY  = node->upper_y[i] * rdir.y - org_rdir.y;
	      const float nearZ = node->lower_z[i] * rdir.z - org_rdir.z;
	      const float farZ  = node->upper_z[i] * rdir.z - org_rdir.z;
	      const float tNearX = min(nearX,farX);
	      const float tFarX  = max(nearX,farX);
	      const float tNearY = min(nearY,farY);
	      const float tFarY  = max(nearY,farY);
	      const float tNearZ = min(nearZ,farZ);
	      const float tFarZ  = max(nearZ,farZ);
          
	      const float tNear = max(tNearX,tNearY,tNearZ,ray.tnear);
	      const float tFar  = min(tFarX ,tFarY ,tFarZ ,ray.tfar);
	      if (tNear <= tFar)
		{
		  stackPtr->ptr  = node->child(i);
		  stackPtr->dist = tNear;
		  stackPtr++;
		  pushed++;
		}
	    }

	  if (pushed == 0) 
	    {
	      goto pop;
	    }
	  else if (pushed == 1)
	    {
	      cur = (NodeRef) stackPtr[-1].ptr;
	      stackPtr--;
	      continue;
	    }
	  else if (pushed == 2)
	    {
	      sort(stackPtr[-1],stackPtr[-2]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;
	      continue;	      
	    }
	  else if (pushed == 3)
	    {
	      sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;
	      continue;	      
	    }
	  else
	    {
	      sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;	      
	    }
	  
        }
        
        /*! this is a leaf node */
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Triangle1* tri = (Triangle1*) cur.leaf(triPtr,num);
        //TriangleIntersector::intersect(ray,tri,num,bvh->geometry);
	for (size_t i=0;i<num;i++)
	  intersect_vec3f(ray,tri[i],bvh->geometry);
      }
      AVX_ZERO_UPPER();
    }
    
    template<typename TriangleIntersector>
    void BVH4iIntersector1Scalar<TriangleIntersector>::occluded(const BVH4i* bvh, Ray& ray)
    {
      /*! stack state */
      StackItem stack[1+3*BVH4i::maxDepth];  //!< stack of nodes 
      StackItem *stackPtr = stack+1;        //!< current stack pointer
      stack[0].ptr  = bvh->root;
      stack[0].dist = neg_inf;
      
      
      /*! load the ray into SIMD registers */
      const Vec3f rdir = rcp_safe(ray.dir);
      const Vec3f org_rdir = ray.org*rdir;
      
      const void* nodePtr = bvh->nodePtr();
      const void* triPtr  = bvh->triPtr();
      
      /* pop loop */
      while (true) pop:
      {
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        NodeRef cur = NodeRef(stackPtr->ptr);
        
        /* downtraversal loop */
        while (true)
        {
          /*! stop if we found a leaf */
          if (unlikely(cur.isLeaf())) break;
          STAT3(shadow.trav_nodes,1,1,1);
          
          /*! single ray intersection with 4 boxes */
          const Node* node = cur.node(nodePtr);

	  size_t pushed = 0;
	  for (size_t i=0;i<4;i++)
	    {
	      const float nearX = node->lower_x[i] * rdir.x - org_rdir.x;
	      const float farX  = node->upper_x[i] * rdir.x - org_rdir.x;
	      const float nearY = node->lower_y[i] * rdir.y - org_rdir.y;
	      const float farY  = node->upper_y[i] * rdir.y - org_rdir.y;
	      const float nearZ = node->lower_z[i] * rdir.z - org_rdir.z;
	      const float farZ  = node->upper_z[i] * rdir.z - org_rdir.z;
	      const float tNearX = min(nearX,farX);
	      const float tFarX  = max(nearX,farX);
	      const float tNearY = min(nearY,farY);
	      const float tFarY  = max(nearY,farY);
	      const float tNearZ = min(nearZ,farZ);
	      const float tFarZ  = max(nearZ,farZ);
          
	      const float tNear = max(tNearX,tNearY,tNearZ,ray.tnear);
	      const float tFar  = min(tFarX ,tFarY ,tFarZ ,ray.tfar);

	      if (tNear <= tFar)
		{
		  stackPtr->ptr  = node->child(i);
		  stackPtr->dist = tNear;
		  stackPtr++;
		  pushed++;
		}
	      
	    }

	  if (pushed == 0) 
	    {
	      goto pop;
	    }
	  else if (pushed == 1)
	    {
	      cur = (NodeRef) stackPtr[-1].ptr;
	      stackPtr--;
	      continue;
	    }
	  else if (pushed == 2)
	    {
	      sort(stackPtr[-1],stackPtr[-2]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;
	      continue;	      
	    }
	  else if (pushed == 3)
	    {
	      sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;
	      continue;	      
	    }
	  else
	    {
	      sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
	      cur = (NodeRef) stackPtr[-1].ptr; 
	      stackPtr--;	      
	    }
          
        }
        
        /*! this is a leaf node */
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle1* tri = (Triangle1*) cur.leaf(triPtr,num);
	for (size_t i=0;i<num;i++)
	  if (occluded_vec3f(ray,tri[i],bvh->geometry)) 
	    {
	      ray.geomID = 0;
	      break;
	    }
	if (ray.geomID == 0) break;	
      }
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR1(BVH4iTriangle1Intersector1ScalarMoeller,BVH4iIntersector1Scalar<Triangle1Intersector1MoellerTrumbore>);
    DEFINE_INTERSECTOR1(BVH4iVirtualIntersector1Scalar,BVH4iIntersector1Scalar<VirtualAccelIntersector1>);
  }
}
