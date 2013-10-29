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

#ifndef __EMBREE_ACCEL_TRIANGLE1I_H__
#define __EMBREE_ACCEL_TRIANGLE1I_H__

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Single triangle from an indexed face set. */
  struct Triangle1i
  {
    /*! Name of triangle */
    static const char* name() { return "triangle1i"; }
    
    /*! block size */
    static const size_t blockSize = 1;
    static const size_t logBlockSize = 0;

    /*! Do we need the vertex array. */
    static const bool needVertices = true;

    /*! Cost of ray/triangle intersection. */
    static const int intCost = 1;

    /*! virtual interface to query information about the triangle type */
    static const struct Type : public TriangleType
    {
      Type () : TriangleType("triangle1i",sizeof(Triangle1i),1,true,1) {}

      size_t blocks(size_t x) const {
        return x;
      }
      
      size_t size(const char* This) const {
        return ((Triangle1i*)This)->size();
      }
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle1i*)This)->area((const TriangleMesh*)geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle1i*)This)->pack(prims,(const TriangleMesh*)geom);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle1i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1i (const int& v0, const int& v1, const int& v2, const int& id0, const int& id1)
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) {}

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { return 1; }

     /*! Computes the area of the triangle. */
    __forceinline float area(const TriangleMesh* geom) {
      const Vector3f e1 = geom->vertex(v0)-geom->vertex(v1);
      const Vector3f e2 = geom->vertex(v2)-geom->vertex(v0);
      const Vector3f Ng = cross(e1,e2);
      return length(Ng);
    }
    
    /*! Packs triangle taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const TriangleMesh* geom)
    {
      const PrimRef& prim = *prims; prims++;
      const TriangleMesh::Triangle& tri = geom->triangle(prim.id());
      new (this) Triangle1i(tri.v0,tri.v1,tri.v2,tri.id0 & 0x7FFFFFFF,tri.id1);
    }

  public:
    int32 v0;      //!< Index of 1st vertex.
    int32 v1;      //!< Index of 2nd vertex.
    int32 v2;      //!< Index of 3rd vertex.
    int32 id0;     //!< 1st user ID.
    int32 id1;     //!< 2nd user ID.
  };

  /*! Intersector1 for triangle1i */
  template<typename Algorithm>
  struct Triangle1iIntersector1
  {
    typedef Triangle1i Triangle;

    static __forceinline void intersect(Ray& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      Algorithm::intersect(ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline bool occluded(const Ray& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      return Algorithm::occluded(ray,p0,p1,p2);
    }

    template<typename Ray>
      static __forceinline void intersect(const size_t k, Ray& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      Algorithm::intersect(k,ray,p0,p1,p2,tri.id0,tri.id1);
    }
    
    template<typename Ray>
      static __forceinline bool occluded(const size_t k, const Ray& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      return Algorithm::occluded(k,ray,p0,p1,p2);
    }
  };

#if defined(__SSE__)

  /*! Intersector4 for triangle1i */
  template<typename Algorithm>
  struct Triangle1iIntersector4
  {
    typedef Triangle1i Triangle;
  
    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),4);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline sseb occluded(const sseb& valid, const Ray4& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),4);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      return Algorithm::occluded(valid,ray,p0,p1,p2);
    }
  };

#endif

#if defined(__AVX__)

  /*! Intersector8 for triangle1i */
  template<typename Algorithm>
  struct Triangle1iIntersector8
  {
    typedef Triangle1i Triangle;
    
    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),8);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),8);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      return Algorithm::occluded(valid,ray,p0,p1,p2);
    }
  };

#endif

#if defined(__MIC__)

  /*! Intersector16 for triangle1i */
  template<typename Algorithm>
  struct Triangle1iIntersector16
  {
    typedef Triangle1i Triangle;
    
    static __forceinline void intersect(const mic_m& valid, Ray16& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),16);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline mic_m occluded(const mic_m& valid, const Ray16& ray, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),16);
      const Vector3f& p0 = vertices[tri.v0];
      const Vector3f& p1 = vertices[tri.v1];
      const Vector3f& p2 = vertices[tri.v2];
      return Algorithm::occluded(valid,ray,p0,p1,p2);
    }
  };

#endif
}

#endif


