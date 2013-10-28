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

#include "sys/platform.h"

#if defined(__AVX__)

#include "bvh4_intersector1_avx.h"
#include "../common/stack_item.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector, bool swapX, bool swapY, bool swapZ>
  __forceinline void intersectT(const BVH4* bvh, Ray& ray)
  {
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef StackItemT<size_t> StackItem;
    typedef typename BVH4::NodeRef NodeRef;
    typedef typename BVH4::Node Node;

    /*! stack state */
    StackItem stack[1+3*BVH4::maxDepth];  //!< stack of nodes 
    StackItem* stackPtr = stack+1;        //!< current stack pointer
    stack[0].ptr  = bvh->root;
    stack[0].dist = neg_inf;

    /*! load the ray into SIMD registers */
    const avxf pos_neg = avxf(ssef(+0.0f),ssef(-0.0f));
    const avxf neg_pos = avxf(ssef(-0.0f),ssef(+0.0f));
    const avxf flipSignX = swapX ? neg_pos : pos_neg;
    const avxf flipSignY = swapY ? neg_pos : pos_neg;
    const avxf flipSignZ = swapZ ? neg_pos : pos_neg;
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const avx3f rdir(ray_rdir.x^flipSignX,ray_rdir.y^flipSignY,ray_rdir.z^flipSignZ);
    const avx3f org_rdir(avx3f(ray.org.x,ray.org.y,ray.org.z)*rdir);
    avxf rayNearFar(ssef(ray.tnear),-ssef(ray.tfar));

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

#if defined (__AVX2__) || defined(__MIC__)
        const avxf tLowerUpperX = msub(avxf::load(&node->lower_x), rdir.x, org_rdir.x);
        const avxf tLowerUpperY = msub(avxf::load(&node->lower_y), rdir.y, org_rdir.y);
        const avxf tLowerUpperZ = msub(avxf::load(&node->lower_z), rdir.z, org_rdir.z);
#else
        const avxf tLowerUpperX = (norg.x + avxf::load(&node->lower_x)) * rdir.x;
        const avxf tLowerUpperY = (norg.y + avxf::load(&node->lower_y)) * rdir.y;
        const avxf tLowerUpperZ = (norg.z + avxf::load(&node->lower_z)) * rdir.z;
#endif
        const avxf tNearFarX = swapX ? shuffle<1,0>(tLowerUpperX) : tLowerUpperX;
        const avxf tNearFarY = swapY ? shuffle<1,0>(tLowerUpperY) : tLowerUpperY;
        const avxf tNearFarZ = swapZ ? shuffle<1,0>(tLowerUpperZ) : tLowerUpperZ;
        const avxf tNearFar = max(tNearFarX,tNearFarY,tNearFarZ,rayNearFar);
        const ssef tNear = extract<0>(tNearFar);
        const ssef tFar  = extract<1>(tNearFar);
        size_t mask = movemask(-tNear >= tFar);
                
        /*! if no child is hit, pop next node */
        if (unlikely(mask == 0))
          goto pop;

        /*! one child is hit, continue with that child */
        size_t r = __bsf(mask); mask = __btc(mask,r);
        if (likely(mask == 0)) {
          cur = node->child(r);
          continue;
        }

        /*! two children are hit, push far child, and continue with closer child */
        NodeRef c0 = node->child(r); const float d0 = tNear[r];
        r = __bsf(mask); mask = __btc(mask,r);
        NodeRef c1 = node->child(r); const float d1 = tNear[r];
        if (likely(mask == 0)) {
          if (d0 < d1) { stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++; cur = c0; continue; }
          else         { stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++; cur = c1; continue; }
        }

        /*! Here starts the slow path for 3 or 4 hit children. We push
         *  all nodes onto the stack to sort them there. */
        stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++;
        stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++;

        /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
        r = __bsf(mask); mask = __btc(mask,r);
        NodeRef c = node->child(r); float d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (likely(mask == 0)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
          continue;
        }

        /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
        cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
      }

      /*! this is a leaf node */
      STAT3(normal.trav_leaves,1,1,1);
      size_t num; Triangle* tri = (Triangle*) cur.leaf(triPtr,num);
      for (size_t i=0; i<num; i++)
        TriangleIntersector::intersect(ray,tri[i],bvh->vertices);
      
      rayNearFar = insert<1>(rayNearFar,-ssef(ray.tfar));
    }
  }

  template<typename TriangleIntersector>
  void BVH4Intersector1AVX<TriangleIntersector>::intersect(const BVH4Intersector1AVX* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(normal.travs,1,1,1);

    const BVH4* bvh = This->bvh;
    int swapX = ray.dir.x < 0.0f;
    int swapY = ray.dir.y < 0.0f;
    int swapZ = ray.dir.z < 0.0f;
    int swap = 4*swapX+2*swapY+swapZ;
    
    switch (swap) {
    case 0: intersectT<TriangleIntersector,false,false,false>(bvh,ray); break;
    case 1: intersectT<TriangleIntersector,false,false,true >(bvh,ray); break;
    case 2: intersectT<TriangleIntersector,false,true ,false>(bvh,ray); break;
    case 3: intersectT<TriangleIntersector,false,true ,true >(bvh,ray); break;
    case 4: intersectT<TriangleIntersector,true ,false,false>(bvh,ray); break;
    case 5: intersectT<TriangleIntersector,true ,false,true >(bvh,ray); break;
    case 6: intersectT<TriangleIntersector,true ,true ,false>(bvh,ray); break;
    case 7: intersectT<TriangleIntersector,true ,true ,true >(bvh,ray); break;
    }

    AVX_ZERO_UPPER();
  }

  template<typename TriangleIntersector, bool swapX, bool swapY, bool swapZ>
  __forceinline bool occludedT(const BVH4* bvh, Ray& ray)
  {
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef StackItemT<size_t> StackItem;
    typedef typename BVH4::NodeRef NodeRef;
    typedef typename BVH4::Node Node;

    /*! stack state */
    NodeRef stack[1+3*BVH4::maxDepth];  //!< stack of nodes that still need to get traversed
    NodeRef* stackPtr = stack+1;        //!< current stack pointer
    stack[0]  = bvh->root;
  
    /*! load the ray into SIMD registers */
    const avxf pos_neg = avxf(ssef(+0.0f),ssef(-0.0f));
    const avxf neg_pos = avxf(ssef(-0.0f),ssef(+0.0f));
    const avxf flipSignX = swapX ? neg_pos : pos_neg;
    const avxf flipSignY = swapY ? neg_pos : pos_neg;
    const avxf flipSignZ = swapZ ? neg_pos : pos_neg;
    const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const avx3f rdir(ray_rdir.x^flipSignX,ray_rdir.y^flipSignY,ray_rdir.z^flipSignZ);
    const avx3f org_rdir(avx3f(ray.org.x,ray.org.y,ray.org.z)*rdir);
    const avxf rayNearFar(ssef(ray.tnear),-ssef(ray.tfar));

    const void* nodePtr = bvh->nodePtr();
    const void* triPtr  = bvh->triPtr();
    
    /* pop loop */
    while (true) pop:
    {
      /*! pop next node */
      if (unlikely(stackPtr == stack)) break;
      stackPtr--;
      NodeRef cur = (NodeRef) *stackPtr;
      
      /* downtraversal loop */
      while (true)
      {
        /*! stop if we found a leaf */
        if (unlikely(cur.isLeaf())) break;
        STAT3(shadow.trav_nodes,1,1,1);
	
        /*! single ray intersection with 4 boxes */
        const Node* node = cur.node(nodePtr);
       
#if defined (__AVX2__) || defined(__MIC__)
        const avxf tLowerUpperX = msub(avxf::load(&node->lower_x), rdir.x, org_rdir.x);
        const avxf tLowerUpperY = msub(avxf::load(&node->lower_y), rdir.y, org_rdir.y);
        const avxf tLowerUpperZ = msub(avxf::load(&node->lower_z), rdir.z, org_rdir.z);
#else
        const avxf tLowerUpperX = (norg.x + avxf::load(&node->lower_x)) * rdir.x;
        const avxf tLowerUpperY = (norg.y + avxf::load(&node->lower_y)) * rdir.y;
        const avxf tLowerUpperZ = (norg.z + avxf::load(&node->lower_z)) * rdir.z;
#endif
        const avxf tNearFarX = swapX ? shuffle<1,0>(tLowerUpperX) : tLowerUpperX;
        const avxf tNearFarY = swapY ? shuffle<1,0>(tLowerUpperY) : tLowerUpperY;
        const avxf tNearFarZ = swapZ ? shuffle<1,0>(tLowerUpperZ) : tLowerUpperZ;
        const avxf tNearFar = max(tNearFarX,tNearFarY,tNearFarZ,rayNearFar);
        const ssef tNear = extract<0>(tNearFar);
        const ssef tFar  = extract<1>(tNearFar);
        size_t mask = movemask(-tNear >= tFar);
        
        /*! if no child is hit, pop next node */
        if (unlikely(mask == 0))
          goto pop;

        /*! one child is hit, continue with that child */
        size_t r = __bsf(mask); mask = __btc(mask,r);
        if (likely(mask == 0)) {
          cur = node->child(r);
          continue;
        }

        /*! two children are hit, push far child, and continue with closer child */
        NodeRef c0 = node->child(r); const float d0 = tNear[r];
        r = __bsf(mask); mask = __btc(mask,r);
        NodeRef c1 = node->child(r); const float d1 = tNear[r];
        if (likely(mask == 0)) {
          if (d0 < d1) { *stackPtr = c1; stackPtr++; cur = c0; continue; }
          else         { *stackPtr = c0; stackPtr++; cur = c1; continue; }
        }
        *stackPtr = c0; stackPtr++;
        *stackPtr = c1; stackPtr++;

        /*! three children are hit */
        r = __bsf(mask); mask = __btc(mask,r);
        cur = node->child(r); *stackPtr = cur; stackPtr++;
        if (likely(mask == 0)) {
          stackPtr--;
          continue;
        }

        /*! four children are hit */
        cur = node->child(3);
      }

      /*! this is a leaf node */
      STAT3(shadow.trav_leaves,1,1,1);
      size_t num; Triangle* tri = (Triangle*) cur.leaf(triPtr,num);
      for (size_t i=0; i<num; i++) {
        if (TriangleIntersector::occluded(ray,tri[i],bvh->vertices)) {
          AVX_ZERO_UPPER();
          return true;
        }
      }
    }
   
    AVX_ZERO_UPPER();
    return false;
  }


  template<typename TriangleIntersector>
  bool BVH4Intersector1AVX<TriangleIntersector>::occluded(const BVH4Intersector1AVX* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(shadow.travs,1,1,1);

    const BVH4* bvh = This->bvh;
    int swapX = ray.dir.x < 0.0f;
    int swapY = ray.dir.y < 0.0f;
    int swapZ = ray.dir.z < 0.0f;
    int swap = 4*swapX+2*swapY+swapZ;
    
    switch (swap) {
    case 0: return occludedT<TriangleIntersector,false,false,false>(bvh,ray); break;
    case 1: return occludedT<TriangleIntersector,false,false,true >(bvh,ray); break;
    case 2: return occludedT<TriangleIntersector,false,true ,false>(bvh,ray); break;
    case 3: return occludedT<TriangleIntersector,false,true ,true >(bvh,ray); break;
    case 4: return occludedT<TriangleIntersector,true ,false,false>(bvh,ray); break;
    case 5: return occludedT<TriangleIntersector,true ,false,true >(bvh,ray); break;
    case 6: return occludedT<TriangleIntersector,true ,true ,false>(bvh,ray); break;
    case 7: return occludedT<TriangleIntersector,true ,true ,true >(bvh,ray); break;
    default: return false;
    }
  }

  void BVH4Intersector1AVXRegister () 
  {
    TriangleMesh::intersectors1.add("bvh4","triangle1i","avx","moeller" ,false,BVH4Intersector1AVX<Triangle1iIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle1i","avx","pluecker",true ,BVH4Intersector1AVX<Triangle1iIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle1v","avx","moeller" ,false,BVH4Intersector1AVX<Triangle1vIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle1v","avx","pluecker",true ,BVH4Intersector1AVX<Triangle1vIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle4i","avx","moeller" ,false,BVH4Intersector1AVX<Triangle4iIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle4i","avx","pluecker",true ,BVH4Intersector1AVX<Triangle4iIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle4v","avx","moeller" ,false,BVH4Intersector1AVX<Triangle4vIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle4v","avx","pluecker",true ,BVH4Intersector1AVX<Triangle4vIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh4","triangle1" ,"avx","moeller" ,true ,BVH4Intersector1AVX<Triangle1Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.add("bvh4","triangle4" ,"avx","moeller" ,true ,BVH4Intersector1AVX<Triangle4Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.add("bvh4","triangle8" ,"avx","moeller" ,true ,BVH4Intersector1AVX<Triangle8Intersector1MoellerTrumbore>::create);
    
    VirtualScene::intersectors1.add("bvh4","virtual"   ,"avx","virtual" ,true ,BVH4Intersector1AVX<VirtualObjectIntersector1>::create);
  }
}
#endif

