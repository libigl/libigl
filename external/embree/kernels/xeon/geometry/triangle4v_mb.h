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
  /*! Stores the vertices of 4 triangles in struct of array layout. */
  struct Triangle4vMB
  {
    struct Type : public PrimitiveType 
    {
      Type ();
      size_t blocks(size_t x) const;
      size_t size(const char* This) const;
      std::pair<BBox3fa,BBox3fa> update2(char* prim, size_t num, void* geom) const;
    };

    static Type type;

  public:

    /*! Default constructor. */
    __forceinline Triangle4vMB () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4vMB (const sse3f& a0, const sse3f& a1, 
				const sse3f& b0, const sse3f& b1,
				const sse3f& c0, const sse3f& c1, 
				const ssei& geomIDs, const ssei& primIDs, const ssei& mask, const bool last)
      : v0(a0), v1(b0), v2(c0), d0(a1-a0), d1(b1-b0), d2(c1-c0), geomIDs(geomIDs), primIDs(primIDs | (last << 31))
    {
#if defined(RTCORE_RAY_MASK)
      this->mask = mask;
#endif
    }

    /*! Returns if the specified triangle is valid. */
    __forceinline bool valid(const size_t i) const { 
      assert(i<4); 
      return geomIDs[i] != -1; 
    }

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return geomIDs != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return bitscan(~movemask(valid()));
    }

    /*! calculate the bounds of the triangles at t0 */
    __forceinline BBox3fa bounds0() const 
    {
      sse3f lower = min(v0,v1,v2);
      sse3f upper = max(v0,v1,v2);
      const sseb mask = valid();
      lower.x = select(mask,lower.x,ssef(pos_inf));
      lower.y = select(mask,lower.y,ssef(pos_inf));
      lower.z = select(mask,lower.z,ssef(pos_inf));
      upper.x = select(mask,upper.x,ssef(neg_inf));
      upper.y = select(mask,upper.y,ssef(neg_inf));
      upper.z = select(mask,upper.z,ssef(neg_inf));
      return BBox3fa(Vec3fa(reduce_min(lower.x),reduce_min(lower.y),reduce_min(lower.z)),
		     Vec3fa(reduce_max(upper.x),reduce_max(upper.y),reduce_max(upper.z)));
    }

    /*! calculate the bounds of the triangles at t1 */
    __forceinline BBox3fa bounds1() const 
    {
      const sse3f p0 = v0+d0;
      const sse3f p1 = v1+d1;
      const sse3f p2 = v2+d2;
      sse3f lower = min(p0,p1,p2);
      sse3f upper = max(p0,p1,p2);
      const sseb mask = valid();
      lower.x = select(mask,lower.x,ssef(pos_inf));
      lower.y = select(mask,lower.y,ssef(pos_inf));
      lower.z = select(mask,lower.z,ssef(pos_inf));
      upper.x = select(mask,upper.x,ssef(neg_inf));
      upper.y = select(mask,upper.y,ssef(neg_inf));
      upper.z = select(mask,upper.z,ssef(neg_inf));
      return BBox3fa(Vec3fa(reduce_min(lower.x),reduce_min(lower.y),reduce_min(lower.z)),
		     Vec3fa(reduce_max(upper.x),reduce_max(upper.y),reduce_max(upper.z)));
    }

    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return (N+3)/4; }

    /*! checks if this is the last triangle in the list */
    __forceinline int last() const { 
      return primIDs[0] & 0x80000000; 
    }

    /*! returns the geometry IDs */
    template<bool list>
    __forceinline ssei geomID() const { 
      return geomIDs; 
    }
    template<bool list>
    __forceinline int geomID(const size_t i) const { 
      assert(i<4); return geomIDs[i]; 
    }

    /*! returns the primitive IDs */
    template<bool list>
    __forceinline ssei primID() const { 
      if (list) return primIDs & 0x7FFFFFFF; 
      else      return primIDs;
    }
    template<bool list>
    __forceinline int  primID(const size_t i) const { 
      assert(i<4); 
      if (list) return primIDs[i] & 0x7FFFFFFF; 
      else      return primIDs[i];
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f va0 = zero, vb0 = zero, vc0 = zero;
      sse3f va1 = zero, vb1 = zero, vc1 = zero;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
	const PrimRef& prim = *prims;
	const size_t geomID = prim.geomID();
        const size_t primID = prim.primID();
        const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
	const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
        const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
	const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
        const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
	const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        va0.x[i] = a0.x; va0.y[i] = a0.y; va0.z[i] = a0.z;
	va1.x[i] = a1.x; va1.y[i] = a1.y; va1.z[i] = a1.z;
	vb0.x[i] = b0.x; vb0.y[i] = b0.y; vb0.z[i] = b0.z;
	vb1.x[i] = b1.x; vb1.y[i] = b1.y; vb1.z[i] = b1.z;
	vc0.x[i] = c0.x; vc0.y[i] = c0.y; vc0.z[i] = c0.z;
	vc1.x[i] = c1.x; vc1.y[i] = c1.y; vc1.z[i] = c1.z;
      }
      new (this) Triangle4vMB(va0,va1,vb0,vb1,vc0,vc1,vgeomID,vprimID,vmask,list && !prims); // FIXME: store_nt
    }
    
    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& begin, size_t end, Scene* scene, const bool list)
    {
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f va0 = zero, vb0 = zero, vc0 = zero;
      sse3f va1 = zero, vb1 = zero, vc1 = zero;
      
      for (size_t i=0; i<4 && begin<end; i++, begin++)
      {
	const PrimRef& prim = prims[begin];
        const size_t geomID = prim.geomID();
        const size_t primID = prim.primID();
        const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	const Vec3fa& a0 = mesh->vertex(tri.v[0],0);
	const Vec3fa& a1 = mesh->vertex(tri.v[0],1);
        const Vec3fa& b0 = mesh->vertex(tri.v[1],0);
	const Vec3fa& b1 = mesh->vertex(tri.v[1],1);
        const Vec3fa& c0 = mesh->vertex(tri.v[2],0);
	const Vec3fa& c1 = mesh->vertex(tri.v[2],1);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        va0.x[i] = a0.x; va0.y[i] = a0.y; va0.z[i] = a0.z;
	va1.x[i] = a1.x; va1.y[i] = a1.y; va1.z[i] = a1.z;
	vb0.x[i] = b0.x; vb0.y[i] = b0.y; vb0.z[i] = b0.z;
	vb1.x[i] = b1.x; vb1.y[i] = b1.y; vb1.z[i] = b1.z;
	vc0.x[i] = c0.x; vc0.y[i] = c0.y; vc0.z[i] = c0.z;
	vc1.x[i] = c1.x; vc1.y[i] = c1.y; vc1.z[i] = c1.z;
      }
      new (this) Triangle4vMB(va0,va1,vb0,vb1,vc0,vc1,vgeomID,vprimID,vmask,list && begin>=end);
    }
   
  public:
    sse3f v0;      //!< 1st vertex of the triangles.
    sse3f v1;      //!< 2nd vertex of the triangles.
    sse3f v2;      //!< 3rd vertex of the triangles.
    sse3f d0;      //!< difference vector between time steps t0 and t1 for first vertex
    sse3f d1;      //!< difference vector between time steps t0 and t1 for second vertex
    sse3f d2;      //!< difference vector between time steps t0 and t1 for third vertex
    ssei geomIDs;  //!< user geometry ID
    ssei primIDs;  //!< primitive ID
#if defined(RTCORE_RAY_MASK)
    ssei mask;     //!< geometry mask
#endif
  };
}
