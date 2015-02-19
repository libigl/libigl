// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "primitive.h"

namespace embree
{
  struct Triangle1v
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle1v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1v (const Vec3fa& v0, const Vec3fa& v1, const Vec3fa& v2, const unsigned geomID, const unsigned primID, const unsigned mask, const bool last)
      : v0(v0,primID | (last << 31)), v1(v1,geomID), v2(v2,mask) {}

    /*! calculate the bounds of the triangle */
    __forceinline BBox3fa bounds() const {
      return merge(BBox3fa(v0),BBox3fa(v1),BBox3fa(v2));
    }

    /*! access hidden members */
    template<bool list>
    __forceinline unsigned primID() const { 
      if (list) return v0.a & 0x7FFFFFFF; 
      else      return v0.a; 
    }
    template<bool list>
    __forceinline unsigned geomID() const { 
      return v1.a; 
    }
    __forceinline unsigned mask  () const { return v2.a; }
    __forceinline int last() const { 
      return v0.a & 0x80000000; 
    }

    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return N; }

    /*! fill triangle from triangle list */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      const PrimRef& prim = *prims;
      prims++;

      const unsigned last = list && !prims;
      const size_t geomID = prim.geomID();
      const size_t primID = prim.primID();
      const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
      const TriangleMesh::Triangle& tri = mesh->triangle(primID);
      
      const ssef p0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
      const ssef p1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
      const ssef p2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
      
      store4f_nt(&v0,cast(insert<3>(cast(p0),primID | (last << 31))));
      store4f_nt(&v1,cast(insert<3>(cast(p1),geomID)));
      store4f_nt(&v2,cast(insert<3>(cast(p2),mesh->mask)));
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& i, size_t end, Scene* scene, const bool list)
    {
      const PrimRef& prim = prims[i];
      i++;

      const unsigned last = list && i >= end;
      const size_t geomID = prim.geomID();
      const size_t primID = prim.primID();
      const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
      const TriangleMesh::Triangle& tri = mesh->triangle(primID);
      
      const ssef p0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
      const ssef p1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
      const ssef p2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
      
      store4f_nt(&v0,cast(insert<3>(cast(p0),primID | (last << 31))));
      store4f_nt(&v1,cast(insert<3>(cast(p1),geomID)));
      store4f_nt(&v2,cast(insert<3>(cast(p2),mesh->mask)));
    }

  public:
    Vec3fa v0;          //!< first vertex and primitive ID
    Vec3fa v1;          //!< second vertex and geometry ID
    Vec3fa v2;          //!< third vertex and geometry mask
  };

  struct Triangle1vType : public PrimitiveType {
    static Triangle1vType type;
    Triangle1vType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct TriangleMeshTriangle1v : public Triangle1vType 
  {
    static TriangleMeshTriangle1v type;
    BBox3fa update(char* prim, size_t num, void* geom) const;
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
                                const unsigned geomID, const unsigned primID, const unsigned mask, const bool last)
      : v0(a0,primID | (last << 31)), v1(b0,geomID), v2(c0,mask), d0(a1-a0), d1(b1-b0), d2(c1-c0) {}

    /*! access hidden members */
    template<bool list>
    __forceinline unsigned primID() const { 
      if (list) return v0.a & 0x7FFFFFFF; 
      else      return v0.a; 
    }
    template<bool list>
    __forceinline unsigned geomID() const { 
      return v1.a; 
    }
    __forceinline unsigned mask  () const { return v2.a; }
    __forceinline int last() const { 
      return v0.a & 0x80000000; 
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      const PrimRef& prim = *prims; prims++;
      const unsigned geomID = prim.geomID();
      const unsigned primID = prim.primID();
      const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
      const TriangleMesh::Triangle& tri = mesh->triangle(primID);
      const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
      const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
      const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
      const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
      const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
      const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
      new (this) Triangle1vMB(a0,a1,b0,b1,c0,c1,mesh->id,primID,mesh->mask,list && !prims);
    }
    
  public:
    Vec3fa v0;          //!< first vertex at t0 (and primitive ID)
    Vec3fa v1;          //!< second vertex at t0 (and geometry ID)
    Vec3fa v2;          //!< third vertex at t0 (and geometry mask)
    Vec3fa d0;          //!< difference vector between time steps t0 and t1 for first vertex
    Vec3fa d1;          //!< difference vector between time steps t0 and t1 for second vertex
    Vec3fa d2;          //!< difference vector between time steps t0 and t1 for third vertex
  };

  struct Triangle1vMBType : public PrimitiveType {
    static Triangle1vMBType type;
    Triangle1vMBType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
    std::pair<BBox3fa,BBox3fa> update2(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle1vMB : public Triangle1vMBType
  {
    static TriangleMeshTriangle1vMB type;
    std::pair<BBox3fa,BBox3fa> update2(char* prim, size_t num, void* geom) const;
  };
}
