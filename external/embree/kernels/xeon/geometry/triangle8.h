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

#ifndef __EMBREE_ACCEL_TRIANGLE8_H__
#define __EMBREE_ACCEL_TRIANGLE8_H__

#include "primitive.h"

namespace embree
{
#if defined __AVX__

  /*! Precalculated representation for 8 triangles. Stores for each
      triangle a base vertex, two edges, and the geometry normal to
      speed up intersection calculations. */
  struct Triangle8
  {
  public:

    /*! Default constructor. */
    __forceinline Triangle8 () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle8 (const avx3f& v0, const avx3f& v1, const avx3f& v2, const avxi& geomID, const avxi& primID, const avxi& mask)
      : v0(v0), e1(v0-v1), e2(v2-v0), Ng(cross(e1,e2)), geomID(geomID), primID(primID)
    {
#if defined(__USE_RAY_MASK__)
      this->mask = mask;
#endif
    }

    /*! Returns a mask that tells which triangles are valid. */
    __forceinline avxb valid() const { return geomID != avxi(-1); }

    /*! Returns a mask that tells which triangles are invalid. */
    __forceinline avxb invalid() const { return geomID == avxi(-1); }

    /*! Returns the number of stored triangles. */
    __forceinline unsigned int size() const {
      return __bsf(~movemask(valid()));
    }

    /*! calculate the bounds of the triangle */
    __forceinline BBox3f bounds() const 
    {
      avx3f p0 = v0;
      avx3f p1 = v0+e1;
      avx3f p2 = v0+e2;
      avx3f lower = min(p0,p1,p2);
      avx3f upper = max(p0,p1,p2);
      avxb mask = valid();
      lower.x = select(mask,lower.x,avxf(pos_inf));
      lower.y = select(mask,lower.y,avxf(pos_inf));
      lower.z = select(mask,lower.z,avxf(pos_inf));
      upper.x = select(mask,upper.x,avxf(neg_inf));
      upper.y = select(mask,upper.y,avxf(neg_inf));
      upper.z = select(mask,upper.z,avxf(neg_inf));
      return BBox3f(Vec3fa(reduce_min(lower.x),reduce_min(lower.y),reduce_min(lower.z)),
                    Vec3fa(reduce_max(upper.x),reduce_max(upper.y),reduce_max(upper.z)));
    }

    /*! non temporal store */
    __forceinline static void store_nt(Triangle8* dst, const Triangle8& src)
    {
      store8f_nt(&dst->v0.x,src.v0.x);
      store8f_nt(&dst->v0.y,src.v0.y);
      store8f_nt(&dst->v0.z,src.v0.z);
      store8f_nt(&dst->e1.x,src.e1.x);
      store8f_nt(&dst->e1.y,src.e1.y);
      store8f_nt(&dst->e1.z,src.e1.z);
      store8f_nt(&dst->e2.x,src.e2.x);
      store8f_nt(&dst->e2.y,src.e2.y);
      store8f_nt(&dst->e2.z,src.e2.z);
      store8f_nt(&dst->Ng.x,src.Ng.x);
      store8f_nt(&dst->Ng.y,src.Ng.y);
      store8f_nt(&dst->Ng.z,src.Ng.z);
      store8i_nt(&dst->geomID,src.geomID);
      store8i_nt(&dst->primID,src.primID);
#if defined(__USE_RAY_MASK__)
      store8i_nt(&dst->mask,src.mask);
#endif
    }

  public:
    avx3f v0;      //!< Base vertex of the triangles.
    avx3f e1;      //!< 1st edge of the triangles (v0-v1).
    avx3f e2;      //!< 2nd edge of the triangles (v2-v0).
    avx3f Ng;      //!< Geometry normal of the triangles.
    avxi geomID;   //!< user geometry ID
    avxi primID;   //!< primitive ID
#if defined(__USE_RAY_MASK__)
    avxi mask;     //!< geometry mask
#endif

  };
#endif

#if defined (__AVX__)

    __forceinline std::ostream &operator<<(std::ostream &o, const Triangle8 &tri)
    {
      o << "v0    " << tri.v0 << std::endl;
      o << "e1    " << tri.e1 << std::endl;
      o << "e2    " << tri.e2 << std::endl;
      o << "Ng    " << tri.Ng << std::endl;
      o << "geomID" << tri.geomID << std::endl;
      o << "primID" << tri.primID << std::endl;
      return o;
    }
#endif


#if defined(__TARGET_AVX__)
  struct Triangle8Type : public PrimitiveType {
    Triangle8Type ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneTriangle8 : public Triangle8Type
  {
    static SceneTriangle8 type;
    void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    void pack(char* dst, const PrimRef* prims, size_t num, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };

  struct TriangleMeshTriangle8 : public Triangle8Type
  {
    static TriangleMeshTriangle8 type;
    void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const;
    BBox3f update(char* prim, size_t num, void* geom) const;
  };
#endif
}

#endif
