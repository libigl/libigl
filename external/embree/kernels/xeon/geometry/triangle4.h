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

#ifndef __EMBREE_ACCEL_TRIANGLE4_H__
#define __EMBREE_ACCEL_TRIANGLE4_H__

#include "primitive.h"

namespace embree
{
  /*! Precalculated representation for 4 triangles. Stores for each
      triangle a base vertex, two edges, and the geometry normal to
      speed up intersection calculations. */
  struct Triangle4
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle4 () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4 (const sse3f& v0, const sse3f& v1, const sse3f& v2, const ssei& geomID, const ssei& primID, const ssei& mask)
      : v0(v0), e1(v0-v1), e2(v2-v0), Ng(cross(e1,e2)), geomID(geomID), primID(primID)
    {
#if defined(__USE_RAY_MASK__)
      this->mask = mask;
#endif
    }

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return geomID != ssei(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return bitscan(~movemask(valid()));
    }

    /*! calculate the bounds of the triangle */
    __forceinline BBox3f bounds() const 
    {
      sse3f p0 = v0;
      sse3f p1 = v0+e1;
      sse3f p2 = v0+e2;
      sse3f lower = min(p0,p1,p2);
      sse3f upper = max(p0,p1,p2);
      sseb mask = valid();
      lower.x = select(mask,lower.x,ssef(pos_inf));
      lower.y = select(mask,lower.y,ssef(pos_inf));
      lower.z = select(mask,lower.z,ssef(pos_inf));
      upper.x = select(mask,upper.x,ssef(neg_inf));
      upper.y = select(mask,upper.y,ssef(neg_inf));
      upper.z = select(mask,upper.z,ssef(neg_inf));
      return BBox3f(Vec3fa(reduce_min(lower.x),reduce_min(lower.y),reduce_min(lower.z)),
                    Vec3fa(reduce_max(upper.x),reduce_max(upper.y),reduce_max(upper.z)));
    }

    /*! non temporal store */
    __forceinline static void store_nt(Triangle4* dst, const Triangle4& src)
    {
      store4f_nt(&dst->v0.x,src.v0.x);
      store4f_nt(&dst->v0.y,src.v0.y);
      store4f_nt(&dst->v0.z,src.v0.z);
      store4f_nt(&dst->e1.x,src.e1.x);
      store4f_nt(&dst->e1.y,src.e1.y);
      store4f_nt(&dst->e1.z,src.e1.z);
      store4f_nt(&dst->e2.x,src.e2.x);
      store4f_nt(&dst->e2.y,src.e2.y);
      store4f_nt(&dst->e2.z,src.e2.z);
      store4f_nt(&dst->Ng.x,src.Ng.x);
      store4f_nt(&dst->Ng.y,src.Ng.y);
      store4f_nt(&dst->Ng.z,src.Ng.z);
      store4i_nt(&dst->geomID,src.geomID);
      store4i_nt(&dst->primID,src.primID);
#if defined(__USE_RAY_MASK__)
      store4i_nt(&dst->mask,src.mask);
#endif
    }

  public:
    sse3f v0;      //!< Base vertex of the triangles.
    sse3f e1;      //!< 1st edge of the triangles (v0-v1).
    sse3f e2;      //!< 2nd edge of the triangles (v2-v0).
    sse3f Ng;      //!< Geometry normal of the triangles.
    ssei geomID;   //!< user geometry ID
    ssei primID;   //!< primitive ID
#if defined(__USE_RAY_MASK__)
    ssei mask;     //!< geometry mask
#endif
  };

  struct Triangle4Type : public PrimitiveType {
    Triangle4Type ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneTriangle4 : public Triangle4Type
  {
    static SceneTriangle4 type;
    void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle4 : public Triangle4Type
  {
    static TriangleMeshTriangle4 type;
    void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };
}

#endif
