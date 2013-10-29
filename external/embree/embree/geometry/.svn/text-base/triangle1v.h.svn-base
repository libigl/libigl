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

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Triangle representation that stores pre-gathered vertices. */
  struct Triangle1v
  {
    /*! Name of triangle */
    static const char* name() { return "triangle1v"; }

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
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle1v*)This)->area((const TriangleMesh*)geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle1v*)This)->pack(prims,(const TriangleMesh*)geom);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle1v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle1v (const Vector3f& v0, const Vector3f& v1, const Vector3f& v2, const int32& id0, const int32& id1)
      : v0(v0), v1(v1), v2(v2) { this->v0.a = id0; this->v1.a = id1; }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return 1;
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const TriangleMesh* geom) {
      return length(cross(v0-v1,v2-v0));
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
      new (this) Triangle1v(p0,p1,p2,tri.id0 & 0x7FFFFFFF,tri.id1);
    }

  public:
    Vec3fa v0;      //!< First vertex of triangle and 1st user ID
    Vec3fa v1;      //!< Second vertex of triangle and 2nd user ID
    Vec3fa v2;      //!< Third vertex of triangle
  };

  /*! Intersector1 for triangle1v */
  template<typename Algorithm>
  struct Triangle1vIntersector1
  {
    typedef Triangle1v Triangle;

    static __forceinline void intersect(Ray& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      Algorithm::intersect(ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2,tri.v0.a,tri.v1.a);
    }

    static __forceinline bool occluded(const Ray& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      return Algorithm::occluded(ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2);
    }
  };
  
#if defined(__SSE__)

  /*! Intersector4 for triangle1v */
  template<typename Algorithm>
  struct Triangle1vIntersector4
  {
    typedef Triangle1v Triangle;
  
    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),4);
      Algorithm::intersect(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2,tri.v0.a,tri.v1.a);
    }

    static __forceinline sseb occluded(const sseb& valid, const Ray4& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),4);
      return Algorithm::occluded(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2);
    }
  };

#endif

#if defined(__AVX__)

  /*! Intersector8 for triangle1v */
  template<typename Algorithm>
  struct Triangle1vIntersector8
  {
    typedef Triangle1v Triangle;
    
    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),8);
      Algorithm::intersect(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2,tri.v0.a,tri.v1.a);
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),8);
      return Algorithm::occluded(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2);
    }
  };

#endif

#if defined(__MIC__)

  /*! Intersector16 for triangle1v */
  template<typename Algorithm>
  struct Triangle1vIntersector16
  {
    typedef Triangle1v Triangle;
    
    static __forceinline void intersect(const mic_m& valid, Ray16& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),16);
      Algorithm::intersect(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2,tri.v0.a,tri.v1.a);
    }

    static __forceinline mic_m occluded(const mic_m& valid, const Ray16& ray, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),16);
      return Algorithm::occluded(valid,ray,(Vector3f&)tri.v0,(Vector3f&)tri.v1,(Vector3f&)tri.v2);
    }
  };

  
#endif
}

#endif


