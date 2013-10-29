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

#ifndef __EMBREE_ACCEL_TRIANGLE4_H__
#define __EMBREE_ACCEL_TRIANGLE4_H__

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Precalculated representation for 4 triangles. Stores for each
      triangle a base vertex, two edges, and the geometry normal to
      speed up intersection calculations. */
  struct Triangle4
  {
    /*! Name of triangle */
    static const char* name() { return "triangle4"; }

    /*! block size */
    static const size_t blockSize = 4;
    static const size_t logBlockSize = 2;

    /*! Tests if we need the vertex array. */
    static const bool needVertices = false;

    /*! Cost of ray/triangle intersection. */
    static const int intCost = 1;

    /*! virtual interface to query information about the triangle type */
    static const struct Type : public TriangleType
    {
      Type () : TriangleType("triangle4",sizeof(Triangle4),4,false,1) {} 

      size_t blocks(size_t x) const {
        return (x+3)/4;
      }
      
      size_t size(const char* This) const {
        return ((Triangle4*)This)->size();
      }
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle4*)This)->area((const TriangleMesh*)geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle4*)This)->pack(prims,(const TriangleMesh*)geom);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle4 () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4 (const sse3f& v0, const sse3f& v1, const sse3f& v2, const ssei& id0, const ssei& id1)
      : v0(v0), e1(v0-v1), e2(v2-v0), Ng(cross(e1,e2)), id0(id0), id1(id1) {}

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return id0 != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return __bsf(~movemask(valid()));
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const TriangleMesh* geom) {
      return reduce_add(select(valid(),length(Ng),ssef(0.0f)));
    }

    /*! Packs 4 triangles taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const TriangleMesh* geom)
    {
      ssei id0 = -1, id1 = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
        const PrimRef& prim = *prims;
        const TriangleMesh::Triangle& tri = geom->triangle(prim.id());
        const Vec3fa& p0 = geom->vertex(tri.v0);
        const Vec3fa& p1 = geom->vertex(tri.v1);
        const Vec3fa& p2 = geom->vertex(tri.v2);
        id0 [i] = tri.id0 & 0x7FFFFFFF; // no support for motion blur
        id1 [i] = tri.id1;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      new (this) Triangle4(v0,v1,v2,id0,id1);
    }

  public:
    sse3f_m v0;      //!< Base vertex of the triangles.
    sse3f_m e1;      //!< 1st edge of the triangles (v0-v1).
    sse3f_m e2;      //!< 2nd edge of the triangles (v2-v0).
    sse3f_m Ng;      //!< Geometry normal of the triangles.
    ssei_m id0;      //!< 1st user ID.
    ssei_m id1;      //!< 2nd user ID.
  };
}

#endif


