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

#ifndef __EMBREE_ACCEL_TRIANGLE1V_H__
#define __EMBREE_ACCEL_TRIANGLE1V_H__

#include "triangle.h"

namespace embree
{
  /*! Triangle representation that stores pre-gathered vertices. */
  struct Triangle1v
  {
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
      Type () : TriangleType("triangle1v",sizeof(Triangle1v),1,false,1) {}

      size_t blocks(size_t x) const {
        return x;
      }
      
      size_t size(const char* This) const {
        return ((Triangle1v*)This)->size();
      }
      
      float area(const char* This, const Vec3fa* vertices) const {
        return ((Triangle1v*)This)->area(vertices);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const BuildTriangle* triangles, const Vec3fa* vertices) const {
        ((Triangle1v*)This)->pack(prims,triangles,vertices);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle1v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1v (const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const int32& id0, const int32& id1)
      : v0(v0), v1(v1), v2(v2) { this->v0.a = id0; this->v1.a = id1; }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return 1;
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const Vec3fa* vertices) {
      return length(cross(v0-v1,v2-v0));
    }

    /*! Packs triangle taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const BuildTriangle* triangles, const Vec3fa* vertices)
    {
      const PrimRef& prim = *prims; prims++;
      const BuildTriangle& tri = triangles[prim.id()];
      const Vec3f& p0 = vertices[tri.v0];
      const Vec3f& p1 = vertices[tri.v1];
      const Vec3f& p2 = vertices[tri.v2];
      new (this) Triangle1v(p0,p1,p2,tri.id0 & 0x7FFFFFFF,tri.id1);
    }

  public:
    Vec3fa v0;      //!< First vertex of triangle and 1st user ID
    Vec3fa v1;      //!< Second vertex of triangle and 2nd user ID
    Vec3fa v2;      //!< Third vertex of triangle
  };
}

#endif


