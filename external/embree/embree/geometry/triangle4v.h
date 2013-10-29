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

#ifndef __EMBREE_ACCEL_TRIANGLE4V_H__
#define __EMBREE_ACCEL_TRIANGLE4V_H__

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Triangle representation that stores pre-gathered vertices for 4 triangles. */
  struct Triangle4v
  {
    /*! Name of triangle */
    static const char* name() { return "triangle4v"; }

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
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle4v*)This)->area((const TriangleMesh*) geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle4v*)This)->pack(prims,(const TriangleMesh*) geom);
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
    __forceinline float area(const TriangleMesh* geom) {
      return reduce_add(select(valid(),length(cross(v0-v1,v2-v0)),ssef(0.0f)));
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

      new (this) Triangle4v(v0,v1,v2,id0,id1);
    }

  public:
    sse3f_m v0;       //!< First vertices of triangles
    sse3f_m v1;       //!< Second vertices of triangles
    sse3f_m v2;       //!< Third vertices of triangles
    ssei_m  id0;      //!< 1st user ID.
    ssei_m  id1;      //!< 2nd user ID.
  };

  /*! Intersector1 for triangle4v */
  template<typename Algorithm>
  struct Triangle4vIntersector1
  {
    typedef Triangle4v Triangle;

    static __forceinline void intersect(Ray& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      Algorithm::intersect(ray,tri.v0,tri.v1,tri.v2,tri.id0,tri.id1);
    }

    static __forceinline bool occluded(const Ray& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      return Algorithm::occluded(ray,tri.v0,tri.v1,tri.v2);
    }
  };

#if defined(__SSE__)

  /*! Intersector4 for triangle4v */
  template<typename Algorithm>
  struct Triangle4vIntersector4
  {
    typedef Triangle4v Triangle;
  
    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),4);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline sseb occluded(const sseb& valid, const Ray4& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),4);
      sseb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif

#if defined(__AVX__)

  /*! Intersector8 for triangle4v */
  template<typename Algorithm>
  struct Triangle4vIntersector8
  {
    typedef Triangle4v Triangle;
    
    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),8);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),8);
      avxb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif

#if defined(__MIC__)

  /*! Intersector16 for triangle4v */
  template<typename Algorithm>
  struct Triangle4vIntersector16
  {
    typedef Triangle4v Triangle;
    
    static __forceinline void intersect(const mic_m& valid, Ray16& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),16);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline mic_m occluded(const mic_m& valid, const Ray16& ray, const Triangle4v& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),16);
      mic_m occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f p0(tri.v0.x[i],tri.v0.y[i],tri.v0.z[i]);
        const Vector3f p1(tri.v1.x[i],tri.v1.y[i],tri.v1.z[i]);
        const Vector3f p2(tri.v2.x[i],tri.v2.y[i],tri.v2.z[i]);
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif
}

#endif


