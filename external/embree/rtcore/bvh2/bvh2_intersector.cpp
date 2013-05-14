// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "bvh2_intersector.h"
#include "triangle/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH2Intersector<TriangleIntersector>::intersect(const Ray& ray, Hit& hit) const
  {
    AVX_ZERO_UPPER();
    STAT3(normal.travs,1,1,1);

    struct StackItem {
      Base* ptr;   //!< node pointer
      float dist;  //!< distance of node
    };

    /*! stack state */
    StackItem stack[1+BVH2::maxDepth];  //!< stack of nodes that still need to get traversed
    StackItem* stackPtr = stack;        //!< current stack pointer
    Base* cur = bvh->root;              //!< in cur we track the ID of the current node

    /*! precomputed shuffles, to switch lower and upper bounds depending on ray direction */
    const ssei identity = _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1, 0);
    const ssei swap     = _mm_set_epi8( 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9, 8);
    const ssei shuffleX = ray.dir.x >= 0 ? identity : swap;
    const ssei shuffleY = ray.dir.y >= 0 ? identity : swap;
    const ssei shuffleZ = ray.dir.z >= 0 ? identity : swap;

    /*! load the ray into SIMD registers */
    const ssei pn = ssei(0x00000000,0x00000000,0x80000000,0x80000000);
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const sse3f rdir = sse3f(ssef(ray.rdir.x) ^ pn, ssef(ray.rdir.y) ^ pn, ssef(ray.rdir.z) ^ pn);
    ssef nearFar(ray.near, ray.near, -ray.far, -ray.far);
    hit.t = min(hit.t,ray.far);

    while (true)
    {
      /*! downtraversal loop */
      while (likely(cur->isNode()))
      {
        /*! single ray intersection with box of both children. */
        const Node* node = cur->node();
        const ssef tNearFarX = (shuffle8(node->lower_upper_x,shuffleX) + norg.x) * rdir.x;
        const ssef tNearFarY = (shuffle8(node->lower_upper_y,shuffleY) + norg.y) * rdir.y;
        const ssef tNearFarZ = (shuffle8(node->lower_upper_z,shuffleZ) + norg.z) * rdir.z;
        const ssef tNearFar = max(tNearFarX,tNearFarY,tNearFarZ,nearFar) ^ pn;
        const sseb lrhit = tNearFar <= shuffle8(tNearFar,swap);

        /*! if two children hit, push far node onto stack and continue with closer node */
        if (likely(lrhit[0] != 0 && lrhit[1] != 0)) {
          if (likely(tNearFar[0] < tNearFar[1])) { 
            stackPtr->ptr = node->child[1]; 
            stackPtr->dist = tNearFar[1]; 
            cur = node->child[0]; 
            stackPtr++; 
          }
          else { 
            stackPtr->ptr = node->child[0]; 
            stackPtr->dist = tNearFar[0]; 
            cur = node->child[1]; 
            stackPtr++; 
          }
        }

        /*! if one child hit, continue with that child */
        else {
          if      (likely(lrhit[0] != 0)) cur = node->child[0];
          else if (likely(lrhit[1] != 0)) cur = node->child[1];
          else goto pop_node;
        }
      }

      /*! leaf node, intersect all triangles */
      {
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++)
          TriangleIntersector::intersect(ray,hit,tri[i],bvh->vertices);
        nearFar = shuffle<0,1,2,3>(nearFar,-hit.t);
      }

      /*! pop next node from stack */
pop_node:
      if (unlikely(stackPtr == stack)) break;
      --stackPtr;
      cur = stackPtr->ptr;
      if (unlikely(stackPtr->dist > hit.t)) goto pop_node;
    }
    AVX_ZERO_UPPER();
  }

  template<typename TriangleIntersector>
  bool BVH2Intersector<TriangleIntersector>::occluded(const Ray& ray) const
  {
    AVX_ZERO_UPPER();

    /*! stack state */
    Base* stack[1+BVH2::maxDepth];   //!< stack of nodes that still need to get traversed
    Base** stackPtr = stack;         //!< current stack pointer
    Base* cur = bvh->root;           //!< in cur we track the ID of the current node

    /*! precomputed shuffles, to switch lower and upper bounds depending on ray direction */
    const ssei identity = _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1, 0);
    const ssei swap     = _mm_set_epi8( 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9, 8);
    const ssei shuffleX = ray.dir.x >= 0 ? identity : swap;
    const ssei shuffleY = ray.dir.y >= 0 ? identity : swap;
    const ssei shuffleZ = ray.dir.z >= 0 ? identity : swap;

    /*! load the ray into SIMD registers */
    const ssei pn = ssei(0x00000000,0x00000000,0x80000000,0x80000000);
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const sse3f rdir = sse3f(ssef(ray.rdir.x) ^ pn, ssef(ray.rdir.y) ^ pn, ssef(ray.rdir.z) ^ pn);
    ssef nearFar(ray.near, ray.near, -ray.far, -ray.far);

    while (true)
    {
      /*! this is an inner node */
      while (likely(cur->isNode()))
      {
        /*! Single ray intersection with box of both children. See bvh2i.h for node layout. */
        const Node* node = cur->node();
        const ssef tNearFarX = (shuffle8(node->lower_upper_x,shuffleX) + norg.x) * rdir.x;
        const ssef tNearFarY = (shuffle8(node->lower_upper_y,shuffleY) + norg.y) * rdir.y;
        const ssef tNearFarZ = (shuffle8(node->lower_upper_z,shuffleZ) + norg.z) * rdir.z;
        const ssef tNearFar = max(tNearFarX,tNearFarY,tNearFarZ,nearFar) ^ pn;
        const sseb lrhit = tNearFar <= shuffle8(tNearFar,swap);

        /*! if two children hit, push far node onto stack and continue with closer node */
        if (likely(lrhit[0] != 0 && lrhit[1] != 0)) {
          *stackPtr++ = node->child[0]; cur = node->child[1];
        }

        /*! if one child hit, continue with that child */
        else {
          if      (lrhit[0] != 0) cur = node->child[0];
          else if (lrhit[1] != 0) cur = node->child[1];
          else goto pop_node;
        }
      }

      /*! leaf node, intersect all triangles */
      {
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++)
          if (TriangleIntersector::occluded(ray,tri[i],bvh->vertices)) {
            AVX_ZERO_UPPER();
            return true;
          }
      }

      /*! pop next node from stack */
pop_node:
      if (unlikely(stackPtr == stack)) break;
      cur = *(--stackPtr);
    }
    AVX_ZERO_UPPER();
    return false;
  }

  /* explicit template instantiation */
  INSTANTIATE_TEMPLATE_BY_INTERSECTOR(BVH2Intersector);
}
