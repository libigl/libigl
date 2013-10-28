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

#ifndef __EMBREE_ACCEL_TRIANGLE4I_H__
#define __EMBREE_ACCEL_TRIANGLE4I_H__

#include "triangle.h"
#include "triangle_mesh.h"

namespace embree
{
  /*! Stores 4 triangles from an indexed face set. */
  struct Triangle4i
  {
    /*! Name of triangle */
    static const char* name() { return "triangle4i"; }

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
      Type () : TriangleType("triangle4i",sizeof(Triangle4i),4,true,2) {}

      size_t blocks(size_t x) const {
        return (x+3)/4;
      }
      
      size_t size(const char* This) const {
        return ((Triangle4i*)This)->size();
      }
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((Triangle4i*)This)->area((const TriangleMesh*)geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((Triangle4i*)This)->pack(prims,(const TriangleMesh*) geom);
      }
      
      std::pair<BBox3f,BBox3f> bounds(char* This, const RTCGeometry* geom) const { 
        return ((Triangle4i*)This)->bounds((const TriangleMesh*)geom);
      }

    } type;

  public:

    /*! Default constructor. */
    __forceinline Triangle4i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4i (const ssei& v0, const ssei& v1, const ssei& v2, const ssei& id0, const ssei& id1)
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) {}

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return id0 != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { 
      return __bsf(~movemask(valid()));
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const TriangleMesh* geom) 
    {
      float area = 0.0f;
      for (size_t i=0; i<4; i++) 
      {
        if (id0[i] == -1) continue;
        const Vector3f e1 = geom->vertex(v0[i])-geom->vertex(v1[i]);
        const Vector3f e2 = geom->vertex(v2[i])-geom->vertex(v0[i]);
        const Vector3f Ng = cross(e1,e2);
        area += length(Ng);
      }
      return area;
    }

    /*! Computes the bounds of the i'th, triangle. */
    __forceinline std::pair<BBox3f,BBox3f> bounds(size_t i, const TriangleMesh* geom) 
    {
      /*! moving triangle */
      if (id0[i] & 0x80000000) 
      {
        const Vec3fa& p0 = geom->vertex(v0[i]+0);
        const Vec3fa& p1 = geom->vertex(v1[i]+0);
        const Vec3fa& p2 = geom->vertex(v2[i]+0);
        const BBox3f bounds0 = merge(BBox3f(p0),BBox3f(p1),BBox3f(p2));
        const Vec3fa& dp0 = geom->vertex(v0[i]+1);
        const Vec3fa& dp1 = geom->vertex(v1[i]+1);
        const Vec3fa& dp2 = geom->vertex(v2[i]+1);
        const BBox3f bounds1 = merge(BBox3f(p0+dp0),BBox3f(p1+dp1),BBox3f(p2+dp2));
        return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
      }
      /*! static triangle */
      else
      {
        const Vec3fa& p0 = geom->vertex(v0[i]);
        const Vec3fa& p1 = geom->vertex(v1[i]);
        const Vec3fa& p2 = geom->vertex(v2[i]);
        const BBox3f bounds = merge(BBox3f(p0),BBox3f(p1),BBox3f(p2));
        return std::pair<BBox3f,BBox3f>(bounds,bounds);
      }
    }

    /*! Computes the bounds of the triangle. */
    __forceinline std::pair<BBox3f,BBox3f> bounds(const TriangleMesh* geom) 
    {
      std::pair<BBox3f,BBox3f> b(empty,empty);
      for (size_t i=0; i<size(); i++) {
        std::pair<BBox3f,BBox3f> bi = bounds(i,geom);
        b.first .grow(bi.first );
        b.second.grow(bi.second);
      }
      return b;
    }
    
    /*! Packs 4 triangles taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const TriangleMesh* geom)
    {
      id0 = -1; id1 = -1; 
      v0 = 0; v1 = 0; v2 = 0;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
        const PrimRef& prim = *prims;
        const TriangleMesh::Triangle& tri = geom->triangle(prim.id());
        id0[i] = tri.id0; id1[i] = tri.id1;
        v0[i] = tri.v0; v1[i] = tri.v1; v2[i] = tri.v2;
      }
    }

  public:
    ssei_m v0;       //!< Pointers to 1st vertex.
    ssei_m v1;       //!< Pointers to 2nd vertex.
    ssei_m v2;       //!< Pointers to 3rd vertex.
    ssei_m id0;      //!< 1st user ID.
    ssei_m id1;      //!< 2nd user ID.
  };

  /*! Intersector1 for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector1
  {
    typedef Triangle4i Triangle;

    static __forceinline void intersect(Ray& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      const ssef_m* vert = (const ssef_m*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      Algorithm::intersect(ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline bool occluded(const Ray& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      const ssef_m* vert = (const ssef_m*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      return Algorithm::occluded(ray,p0,p1,p2);
    }
  };

  /*! Intersector1MB for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector1MB
  {
    typedef Triangle4i Triangle;

    static __forceinline void intersect(Ray& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      const ssef_m* vert = (const ssef_m*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      
      const sseb motion = (tri.id0 & ssei(0x80000000)) != ssei(zero);
      if (unlikely(any(tri.valid() & motion))) {
        const ssei dv0 = select(motion,tri.v0+1,tri.v0);
        const ssei dv1 = select(motion,tri.v1+1,tri.v1);
        const ssei dv2 = select(motion,tri.v2+1,tri.v2);
        sse3f dp0; transpose(vert[dv0[0]],vert[dv0[1]],vert[dv0[2]],vert[dv0[3]], dp0.x,dp0.y,dp0.z);
        sse3f dp1; transpose(vert[dv1[0]],vert[dv1[1]],vert[dv1[2]],vert[dv1[3]], dp1.x,dp1.y,dp1.z);
        sse3f dp2; transpose(vert[dv2[0]],vert[dv2[1]],vert[dv2[2]],vert[dv2[3]], dp2.x,dp2.y,dp2.z);
        p0 = select(motion,p0+ssef(ray.time)*dp0,p0);
        p1 = select(motion,p1+ssef(ray.time)*dp1,p1);
        p2 = select(motion,p2+ssef(ray.time)*dp2,p2);
      }
      Algorithm::intersect(ray,p0,p1,p2,tri.id0,tri.id1);
    }

    static __forceinline bool occluded(const Ray& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,1,1);
      const ssef_m* vert = (const ssef_m*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      
      const sseb motion = (tri.id0 & ssei(0x80000000)) != ssei(zero);
      if (unlikely(any(tri.valid() & motion))) {
        const ssei dv0 = select(motion,tri.v0+1,tri.v0);
        const ssei dv1 = select(motion,tri.v1+1,tri.v1);
        const ssei dv2 = select(motion,tri.v2+1,tri.v2);
        sse3f dp0; transpose(vert[dv0[0]],vert[dv0[1]],vert[dv0[2]],vert[dv0[3]], dp0.x,dp0.y,dp0.z);
        sse3f dp1; transpose(vert[dv1[0]],vert[dv1[1]],vert[dv1[2]],vert[dv1[3]], dp1.x,dp1.y,dp1.z);
        sse3f dp2; transpose(vert[dv2[0]],vert[dv2[1]],vert[dv2[2]],vert[dv2[3]], dp2.x,dp2.y,dp2.z);
        p0 = select(motion,p0+ssef(ray.time)*dp0,p0);
        p1 = select(motion,p1+ssef(ray.time)*dp1,p1);
        p2 = select(motion,p2+ssef(ray.time)*dp2,p2);
      }
      return Algorithm::occluded(ray,p0,p1,p2);
    }
  };

#if defined(__SSE__)

   /*! Intersector4 for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector4
  {
    typedef Triangle4i Triangle;
  
    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),4);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline sseb occluded(const sseb& valid, const Ray4& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),4);
      sseb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

  /*! Intersector4MB for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector4MB
  {
    typedef Triangle4i Triangle;
  
    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),4);
      
      for (size_t i=0; i<tri.size(); i++)
      {
        const size_t v0 = tri.v0[i];
        const size_t v1 = tri.v1[i];
        const size_t v2 = tri.v2[i];
        sse3f p0 = sse3f(vertices[v0]);
        sse3f p1 = sse3f(vertices[v1]);
        sse3f p2 = sse3f(vertices[v2]);
        if (tri.id0[i] & 0x80000000) {
          p0 += ray.time*sse3f(vertices[v0+1]);
          p1 += ray.time*sse3f(vertices[v1+1]);
          p2 += ray.time*sse3f(vertices[v2+1]);
        }
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline sseb occluded(const sseb& valid, const Ray4& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),4);

      sseb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const size_t v0 = tri.v0[i];
        const size_t v1 = tri.v1[i];
        const size_t v2 = tri.v2[i];
        sse3f p0 = sse3f(vertices[v0]);
        sse3f p1 = sse3f(vertices[v1]);
        sse3f p2 = sse3f(vertices[v2]);
        if (tri.id0[i] & 0x80000000) {
          p0 += ray.time*sse3f(vertices[v0+1]);
          p1 += ray.time*sse3f(vertices[v1+1]);
          p2 += ray.time*sse3f(vertices[v2+1]);
        }
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif

#if defined(__AVX__)

  /*! Intersector8 for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector8
  {
    typedef Triangle4i Triangle;
    
    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),8);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),8);
      avxb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

  /*! Intersector8MB for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector8MB
  {
    typedef Triangle4i Triangle;
    
    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),8);

      for (size_t i=0; i<tri.size(); i++)
      {
        const size_t v0 = tri.v0[i];
        const size_t v1 = tri.v1[i];
        const size_t v2 = tri.v2[i];
        avx3f p0 = avx3f(vertices[v0]);
        avx3f p1 = avx3f(vertices[v1]);
        avx3f p2 = avx3f(vertices[v2]);
        if (tri.id0[i] & 0x80000000) {
          p0 += ray.time*avx3f(vertices[v0+1]);
          p1 += ray.time*avx3f(vertices[v1+1]);
          p2 += ray.time*avx3f(vertices[v2+1]);
        }
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),8);

      avxb occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const size_t v0 = tri.v0[i];
        const size_t v1 = tri.v1[i];
        const size_t v2 = tri.v2[i];
        avx3f p0 = avx3f(vertices[v0]);
        avx3f p1 = avx3f(vertices[v1]);
        avx3f p2 = avx3f(vertices[v2]);
        if (tri.id0[i] & 0x80000000) {
          p0 += ray.time*avx3f(vertices[v0+1]);
          p1 += ray.time*avx3f(vertices[v1+1]);
          p2 += ray.time*avx3f(vertices[v2+1]);
        }
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif

#if defined(__MIC__)

  /*! Intersector16 for triangle4i */
  template<typename Algorithm>
  struct Triangle4iIntersector16
  {
    typedef Triangle4i Triangle;
    
    static __forceinline void intersect(const mic_m& valid, Ray16& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,popcnt(valid),16);
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        Algorithm::intersect(valid,ray,p0,p1,p2,tri.id0[i],tri.id1[i]);
      }
    }

    static __forceinline mic_m occluded(const mic_m& valid, const Ray16& ray, const Triangle4i& tri, const Vec3fa* vertices)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid),16);
      mic_m occlusion = !valid;
      for (size_t i=0; i<tri.size(); i++)
      {
        const Vector3f& p0 = vertices[tri.v0[i]];
        const Vector3f& p1 = vertices[tri.v1[i]];
        const Vector3f& p2 = vertices[tri.v2[i]];
        occlusion |= Algorithm::occluded(valid,ray,p0,p1,p2);
        if (all(occlusion)) return occlusion;
      }
      return occlusion;
    }
  };

#endif
}

#endif


