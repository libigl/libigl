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

#ifndef __EMBREE_ACCEL_TRIANGLE1V_H__
#define __EMBREE_ACCEL_TRIANGLE1V_H__

#include "primitive.h"

namespace embree
{
  struct Triangle1v
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle1v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1v (const Vec3fa& v0, const Vec3fa& v1, const Vec3fa& v2, const unsigned geomID, const unsigned primID, const unsigned mask)
      : v0(v0,primID), v1(v1,geomID), v2(v2,mask) {}

    /*! calculate the bounds of the triangle */
    __forceinline BBox3f bounds() const {
      return merge(BBox3f(v0),BBox3f(v1),BBox3f(v2));
    }

    /*! access hidden members */
    __forceinline unsigned primID() const { return v0.a; }
    __forceinline unsigned geomID() const { return v1.a; }
    __forceinline unsigned mask  () const { return v2.a; }
    
  public:
    Vec3fa v0;          //!< first vertex and primitive ID
    Vec3fa v1;          //!< second vertex and geometry ID
    Vec3fa v2;          //!< third vertex and geometry mask
  };

  struct Triangle1vType : public PrimitiveType {
    Triangle1vType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneTriangle1v : public Triangle1vType 
  {
    static SceneTriangle1v type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle1v : public Triangle1vType 
  {
    static TriangleMeshTriangle1v type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };

  struct Triangle1vMB
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle1vMB () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1vMB (const Vec3fa& a0, const Vec3fa& a1, 
                                const Vec3fa& b0, const Vec3fa& b1,
                                const Vec3fa& c0, const Vec3fa& c1, 
                                const unsigned geomID, const unsigned primID, const unsigned mask)
      : v0(a0,primID), v1(b0,geomID), v2(c0,mask), d0(a1-a0), d1(b1-b0), d2(c1-c0) {}

    /*! access hidden members */
    __forceinline unsigned primID() const { return v0.a; }
    __forceinline unsigned geomID() const { return v1.a; }
    __forceinline unsigned mask  () const { return v2.a; }
    
  public:
    Vec3fa v0;          //!< first vertex at t0 (and primitive ID)
    Vec3fa v1;          //!< second vertex at t0 (and geometry ID)
    Vec3fa v2;          //!< third vertex at t0 (and geometry mask)
    Vec3fa d0;          //!< difference vector between time steps t0 and t1 for first vertex
    Vec3fa d1;          //!< difference vector between time steps t0 and t1 for second vertex
    Vec3fa d2;          //!< difference vector between time steps t0 and t1 for third vertex
  };

  struct Triangle1vMBType : public PrimitiveType {
    Triangle1vMBType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneTriangle1vMB : public Triangle1vMBType
  {
    static SceneTriangle1vMB type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    std::pair<BBox3f,BBox3f> update2(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle1vMB : public Triangle1vMBType
  {
    static TriangleMeshTriangle1vMB type;
    void pack(char* dst, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    std::pair<BBox3f,BBox3f> update2(char* prim, size_t num, void* geom) const;
  };
}

#endif
