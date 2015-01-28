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

#include "primitive.h"

namespace embree
{
  struct Triangle1
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle1 () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1 (const Vec3fa& v0, const Vec3fa& v1, const Vec3fa& v2, const unsigned int geomID, const unsigned int primID, const unsigned int mask)
      : v0(v0,primID), v1(v1,geomID), v2(v2,mask), Ng(cross(v0-v1,v2-v0)) {}

    /*! calculate the bounds of the triangle */
    __forceinline BBox3f bounds() const {
      return merge(BBox3f(v0),BBox3f(v1),BBox3f(v2));
    }

    /*! access hidden members */
    __forceinline unsigned int primID() const { return v0.a; }
    __forceinline unsigned int geomID() const { return v1.a; }
    __forceinline unsigned int mask  () const { return v2.a; }
    
  public:
    Vec3fa v0;          //!< first vertex and primitive ID
    Vec3fa v1;          //!< second vertex and geometry ID
    Vec3fa v2;          //!< third vertex and geometry mask
    Vec3fa Ng;          //!< Geometry normal of the triangles.
  };

  struct Triangle1Type : public PrimitiveType {
    Triangle1Type ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneTriangle1 : public Triangle1Type
  {
    static SceneTriangle1 type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const; 
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle1 : public Triangle1Type
  {
    static TriangleMeshTriangle1 type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };
}

#endif
