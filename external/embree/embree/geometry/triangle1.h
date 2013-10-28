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

#ifndef __EMBREE_ACCEL_TRIANGLE1_H__
#define __EMBREE_ACCEL_TRIANGLE1_H__

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Precalculated representation for individual triangles. Stores
      base vertex, two edges, and the geometry normal to speed up
      intersection calculations. */
  struct Triangle1
  {
    /*! Name of triangle */
    static const char* name() { return "triangle1"; }

    /*! block size */
    static const size_t blockSize = 1;
    static const size_t logBlockSize = 0;

    /*! Determines if we need the vertex array. */
    static const bool needVertices = false;

    /*! Cost of ray/triangle intersection. */
    static const int intCost = 1;

    /*! virtual interface to query information about the triangle type */
    static const struct Type : public TriangleType
    {
      Type () : TriangleType("triangle1",sizeof(Triangle1),1,false,1) {}

      size_t blocks(size_t x) const {
        return x;
      }
      
      size_t size(const char* This) const {
        return ((Triangle1*)This)->size();
      }
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle1*)This)->area((const TriangleMesh*)geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle1*)This)->pack(prims,(const TriangleMesh*)geom);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle1 () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1 (const Vector3f& v0, const Vector3f& v1, const Vector3f& v2, const int32& id0, const int32& id1)
      : v0(v0), e1(v0-v1), e2(v2-v0), Ng(cross(e1,e2)) { e1.a = id0; e2.a = id1; }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return 1;
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const TriangleMesh* geom) {
      return length(Ng);
    }

    /*! Packs triangle taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const TriangleMesh* geom)
    {
      const PrimRef& prim = *prims; prims++;
      const TriangleMesh::Triangle& tri = geom->triangle(prim.id());
      const Vec3fa& p0 = geom->vertex(tri.v0);
      const Vec3fa& p1 = geom->vertex(tri.v1);
      const Vec3fa& p2 = geom->vertex(tri.v2);
      new (this) Triangle1(p0,p1,p2,tri.id0 & 0x7FFFFFFF,tri.id1);
    }

  public:
    Vec3fa v0;      //!< Base vertex of the triangles.
    Vec3fa e1;     //!< 1st edge of the triangles (v0-v1) and 1st user ID
    Vec3fa e2;     //!< 2nd edge of the triangles (v2-v0) and 2nd user ID
    Vec3fa Ng;      //!< Geometry normal of the triangles.
  };
}

#endif
