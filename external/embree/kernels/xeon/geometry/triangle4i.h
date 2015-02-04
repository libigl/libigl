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
  /*! Stores 4 triangles from an indexed face set. */
  struct Triangle4i
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle4i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4i (Vec3f* base[4], const ssei& v1, const ssei& v2, const ssei& geomIDs, const ssei& primIDs, const bool last)
      : v1(v1), v2(v2), geomIDs(geomIDs | (last << 31)), primIDs(primIDs) // FIXME: use primID for list end tag
      {
        v0[0] = base[0];
        v0[1] = base[1];
        v0[2] = base[2];
        v0[3] = base[3];
      }

    /*! Returns if the specified triangle is valid. */
    __forceinline bool valid(const size_t i) const { 
      assert(i<4); 
      return geomIDs[i] != -1; 
    }

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return primIDs != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { 
      return __bsf(~movemask(valid()));
    }

    /*! calculate the bounds of the triangles */
    __forceinline const BBox3fa bounds() const 
    {
      BBox3fa bounds = empty;
      for (size_t i=0; i<4 && geomIDs[i] != -1; i++)
      {
	const int* base = (int*) v0[i];
	const Vec3fa p0 = *(Vec3fa*)(base);
	const Vec3fa p1 = *(Vec3fa*)(base+v1[i]);
	const Vec3fa p2 = *(Vec3fa*)(base+v2[i]);
	bounds.extend(p0);
	bounds.extend(p1);
	bounds.extend(p2);
      }
      return bounds;
    }

    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return (N+3)/4; }

    /*! checks if this is the last triangle in the list */
    __forceinline int last() const { 
      return geomIDs[0] & 0x80000000; 
    }
    
    /*! returns the geometry IDs */
    template<bool list>
    __forceinline ssei geomID() const { 
      if (list) return geomIDs & 0x7FFFFFFF; 
      else      return geomIDs;
    }
    template<bool list>
    __forceinline int geomID(const size_t i) const { 
      assert(i<4); 
      if (list) return geomIDs[i] & 0x7FFFFFFF; 
      else      return geomIDs[i];
    }

    /*! returns the primitive IDs */
    template<bool list>
    __forceinline ssei primID() const { 
      return primIDs; 
    }
    template<bool list>
    __forceinline int primID(const size_t i) const { 
      assert(i<4); return primIDs[i]; 
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      ssei geomID = -1, primID = -1;
      Vec3f* v0[4] = { NULL, NULL, NULL, NULL };
      ssei v1 = zero, v2 = zero;
      PrimRef& prim = *prims;
      
      for (size_t i=0; i<4; i++)
      {
	const TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
	const TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
	if (prims) {
	  geomID[i] = prim.geomID();
	  primID[i] = prim.primID();
	  v0[i] = (Vec3f*) mesh->vertexPtr(tri.v[0]); 
	  v1[i] = (int*)   mesh->vertexPtr(tri.v[1])-(int*)v0[i]; 
	  v2[i] = (int*)   mesh->vertexPtr(tri.v[2])-(int*)v0[i]; 
	  prims++;
	} else {
	  assert(i);
	  geomID[i] = -1;
	  primID[i] = -1;
	  v0[i] = v0[i-1];
	  v1[i] = 0; 
	  v2[i] = 0;
	}
	if (prims) prim = *prims; 
      }
      
      new (this) Triangle4i(v0,v1,v2,geomID,primID,list && !prims); // FIXME: use non temporal store
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& begin, size_t end, Scene* scene, const bool list)
    {
      ssei geomID = -1, primID = -1;
      Vec3f* v0[4] = { NULL, NULL, NULL, NULL };
      ssei v1 = zero, v2 = zero;
      const PrimRef* prim = &prims[begin];
      
      for (size_t i=0; i<4; i++)
      {
	const TriangleMesh* mesh = scene->getTriangleMesh(prim->geomID());
	const TriangleMesh::Triangle& tri = mesh->triangle(prim->primID());
	if (begin<end) {
	  geomID[i] = prim->geomID();
	  primID[i] = prim->primID();
	  v0[i] = (Vec3f*) mesh->vertexPtr(tri.v[0]); 
	  v1[i] = (int*)   mesh->vertexPtr(tri.v[1])-(int*)v0[i]; 
	  v2[i] = (int*)   mesh->vertexPtr(tri.v[2])-(int*)v0[i]; 
	  begin++;
	} else {
	  assert(i);
	  geomID[i] = -1;
	  primID[i] = -1;
	  v0[i] = v0[i-1];
	  v1[i] = 0; 
	  v2[i] = 0;
	}
	if (begin<end) prim = &prims[begin];
      }
      
      new (this) Triangle4i(v0,v1,v2,geomID,primID,list && begin>=end); // FIXME: use non temporal store
    }

  public:
    const Vec3f* v0[4];  //!< Pointer to 1st vertex.
    ssei v1;             //!< Offset to 2nd vertex.
    ssei v2;             //!< Offset to 3rd vertex.
    ssei geomIDs;        //!< ID of mesh.
    ssei primIDs;        //!< ID of primitive inside mesh.
  };

  /*! virtual interface to query information about the triangle type */
  struct Triangle4iType : public PrimitiveType
  {
    static Triangle4iType type;

    Triangle4iType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct TriangleMeshTriangle4i : public Triangle4iType
  {
    static TriangleMeshTriangle4i type;
    BBox3fa update(char* prim, size_t num, void* geom) const;
  };
}
