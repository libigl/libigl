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

#ifndef __EMBREE_ACCEL_TRIANGLE4I_H__
#define __EMBREE_ACCEL_TRIANGLE4I_H__

#include "triangle.h"

namespace embree
{
  /*! Stores 4 triangles from an indexed face set. */
  struct Triangle4i
  {
    /*! block size */
    static const size_t blockSize = 4;
    static const size_t logBlockSize = 2;

    /*! Tests if we need the vertex array. */
    static const bool needVertices = true;

    /*! Cost of ray/triangle intersection. */
    static const int intCost = 2;

    /*! virtual interface to query information about the triangle type */
    static const struct Type : public TriangleType
    {
      Type () : TriangleType("triangle4i",sizeof(Triangle4i),4,true,2) {}

      size_t blocks(size_t x) const {
        return (x+3)/4;
      }
      
      size_t size(const char* This) const {
        return ((Triangle4i*)This)->size();
      }
      
      float area(const char* This, const Vec3fa* vertices) const {
        return ((Triangle4i*)This)->area(vertices);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const BuildTriangle* triangles, const Vec3fa* vertices) const {
        ((Triangle4i*)This)->pack(prims,triangles,vertices);
      }
      
      std::pair<BBox3f,BBox3f> bounds(char* This, const Vec3fa* vertices) const { 
        return ((Triangle4i*)This)->bounds(vertices);
      }

    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle4i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4i (const ssei& v0, const ssei& v1, const ssei& v2, const ssei& id0, const ssei& id1)
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) {}

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return id0 != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { 
      return __bsf(~movemask(valid()));
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const Vec3fa* vertices) 
    {
      float area = 0.0f;
      for (size_t i=0; i<4; i++) 
      {
        if (id0[i] == -1) continue;
        const Vec3f e1 = vertices[v0[i]]-vertices[v1[i]];
        const Vec3f e2 = vertices[v2[i]]-vertices[v0[i]];
        const Vec3f Ng = cross(e1,e2);
        area += length(Ng);
      }
      return area;
    }

    /*! Computes the bounds of the i'th, triangle. */
    __forceinline std::pair<BBox3f,BBox3f> bounds(size_t i, const Vec3fa* vertices) 
    {
      /*! moving triangle */
      if (id0[i] & 0x80000000) 
      {
        const Vec3f& p0 = vertices[v0[i]+0];
        const Vec3f& p1 = vertices[v1[i]+0];
        const Vec3f& p2 = vertices[v2[i]+0];
        const BBox3f bounds0 = merge(BBox3f(p0),BBox3f(p1),BBox3f(p2));
        const Vec3f& dp0 = vertices[v0[i]+1];
        const Vec3f& dp1 = vertices[v1[i]+1];
        const Vec3f& dp2 = vertices[v2[i]+1];
        const BBox3f bounds1 = merge(BBox3f(p0+dp0),BBox3f(p1+dp1),BBox3f(p2+dp2));
        return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
      }
      /*! static triangle */
      else
      {
        const Vec3f& p0 = vertices[v0[i]];
        const Vec3f& p1 = vertices[v1[i]];
        const Vec3f& p2 = vertices[v2[i]];
        const BBox3f bounds = merge(BBox3f(p0),BBox3f(p1),BBox3f(p2));
        return std::pair<BBox3f,BBox3f>(bounds,bounds);
      }
    }

    /*! Computes the bounds of the triangle. */
    __forceinline std::pair<BBox3f,BBox3f> bounds(const Vec3fa* vertices) 
    {
      std::pair<BBox3f,BBox3f> b(empty,empty);
      for (size_t i=0; i<size(); i++) {
        std::pair<BBox3f,BBox3f> bi = bounds(i,vertices);
        b.first .grow(bi.first );
        b.second.grow(bi.second);
      }
      return b;
    }
    
    /*! Packs 4 triangles taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const BuildTriangle* triangles, const Vec3fa* vertices)
    {
      id0 = -1; id1 = -1; 
      v0 = 0; v1 = 0; v2 = 0;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
        const PrimRef& prim = *prims;
        const BuildTriangle& tri = triangles[prim.id()];
        id0[i] = tri.id0; id1[i] = tri.id1;
        v0[i] = tri.v0; v1[i] = tri.v1; v2[i] = tri.v2;
      }
    }

  public:
    ssei v0;       //!< Pointers to 1st vertex.
    ssei v1;       //!< Pointers to 2nd vertex.
    ssei v2;       //!< Pointers to 3rd vertex.
    ssei id0;      //!< 1st user ID.
    ssei id1;      //!< 2nd user ID.
  };
}

#endif


