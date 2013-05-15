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

#ifndef __EMBREE_ACCEL_TRIANGLE4V_H__
#define __EMBREE_ACCEL_TRIANGLE4V_H__

#include "triangle.h"

namespace embree
{
  /*! Triangle representation that stores pre-gathered vertices for 4 triangles. */
  struct Triangle4v
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
      Type () : TriangleType("triangle4v",sizeof(Triangle4v),4,false,2) {}

      size_t blocks(size_t x) const {
        return (x+3)/4;
      }
      
      size_t size(const char* This) const {
        return ((Triangle4v*)This)->size();
      }
      
      float area(const char* This, const Vec3fa* vertices) const {
        return ((Triangle4v*)This)->area(vertices);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const BuildTriangle* triangles, const Vec3fa* vertices) const {
        ((Triangle4v*)This)->pack(prims,triangles,vertices);
      }
      
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle4v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4v (const sse3f& v0, const sse3f& v1, const sse3f& v2, const ssei& id0, const ssei& id1)
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) {}
    
    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return id0 != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { 
      return __bsf(~movemask(valid()));
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const Vec3fa* vertices) {
      return reduce_add(select(valid(),length(cross(v0-v1,v2-v0)),ssef(0.0f)));
    }

    /*! Packs 4 triangles taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const BuildTriangle* triangles, const Vec3fa* vertices)
    {
      ssei id0 = -1, id1 = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
        const PrimRef& prim = *prims;
        const BuildTriangle& tri = triangles[prim.id()];
        const Vec3f& p0 = vertices[tri.v0];
        const Vec3f& p1 = vertices[tri.v1];
        const Vec3f& p2 = vertices[tri.v2];
        id0 [i] = tri.id0 & 0x7FFFFFFF; // no support for motion blur
        id1 [i] = tri.id1;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }

      new (this) Triangle4v(v0,v1,v2,id0,id1);
    }

  public:
    sse3f v0;       //!< First vertices of triangles
    sse3f v1;       //!< Second vertices of triangles
    sse3f v2;       //!< Third vertices of triangles
    ssei  id0;      //!< 1st user ID.
    ssei  id1;      //!< 2nd user ID.
  };
}

#endif


