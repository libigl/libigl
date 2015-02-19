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

#include "common/default.h"
#include "primitive.h"
#include "bezier1v.h"

namespace embree
{
#if defined(__SSE__) // FIXME: move to other place
  extern ssef sse_coeff0[4];
  extern ssef sse_coeff1[4];
#endif

#if defined(__AVX__)
  extern avxf coeff0[4];
  extern avxf coeff1[4];
#endif

  struct Bezier1i
  {
  public:

    /*! Default constructor. */
    __forceinline Bezier1i () {}

    /*! Construction from vertices and IDs. */
    __forceinline Bezier1i (const unsigned vertexID, const unsigned int geomID, const unsigned int primID, const bool last)
      : vertexID(vertexID), geom(geomID), prim(primID | (last << 31)) {}

    /*! calculate the bounds of the triangle */
    //__forceinline BBox3fa bounds() const {
    //const BBox3fa b = merge(BBox3fa(p[0]),BBox3fa(p[1]),BBox3fa(p[2]),BBox3fa(p[3]));
    //return enlarge(b,Vec3fa(b.upper.w));
    //}

    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return N; }

    template<bool list>
    __forceinline unsigned int primID() const { 
      if (list) return prim & 0x7FFFFFFF; 
      else      return prim;
    }
    template<bool list>
    __forceinline unsigned int geomID() const { 
      return geom; 
    }
    //__forceinline unsigned int mask  () const { return mask; } // FIXME: not implemented yet
    __forceinline int last  () const { 
      return prim & 0x80000000; 
    }
    
    /*! fill from list */
    __forceinline void fill(atomic_set<PrimRefBlockT<BezierPrim> >::block_iterator_unsafe& iter, Scene* scene, const bool list)
    {
      const BezierPrim& curve = *iter; iter++;
      const unsigned geomID = curve.geomID<0>();
      const unsigned primID = curve.primID<0>();
      const BezierCurves* in = (BezierCurves*) scene->get(geomID);
      //const Vec3fa& p0 = in->vertex(in->curve(primID));
      const unsigned int vertexID = in->curve(primID);
      new (this) Bezier1i(vertexID,geomID,primID,list && !iter);
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& i, size_t end, Scene* scene, const bool list)
    {
      const PrimRef& prim = prims[i];
      i++;
      const size_t geomID = prim.geomID();
      const size_t primID = prim.primID();
      const BezierCurves* curves = scene->getBezierCurves(geomID);
      const size_t vertexID = curves->curve(primID);
      //const Vec3fa& p0 = curves->vertex(vertexID+0);
      new (this) Bezier1i(vertexID,geomID,primID,list && i>=end);
    }

  public:
    //const Vec3fa* p;      //!< pointer to first control point (x,y,z,r)
    unsigned int vertexID; //!< index of start vertex
    unsigned int geom;  //!< geometry ID
    unsigned int prim;  //!< primitive ID
  };

  struct Bezier1iMB
  {
  public:

    /*! Default constructor. */
    __forceinline Bezier1iMB () {}

    /*! Construction from vertices and IDs. */
    __forceinline Bezier1iMB (const unsigned vertexID, const unsigned int geomID, const unsigned int primID, const bool last)
      : vertexID(vertexID), geom(geomID), prim(primID | (last << 31)) {}

    template<bool list>
    __forceinline unsigned int primID() const { 
      if (list) return prim & 0x7FFFFFF; 
      else      return prim;
    }
    template<bool list>
    __forceinline unsigned int geomID() const { 
      return geom; 
    }
    //__forceinline unsigned int mask  () const { return mask; } // FIXME: not implemented yet
    __forceinline int last  () const { 
      return prim & 0x80000000; 
    }

    /*! calculate the bounds of the triangle */
    //__forceinline BBox3fa bounds0() const {
    //const BBox3fa b = merge(BBox3fa(p0[0]),BBox3fa(p0[1]),BBox3fa(p0[2]),BBox3fa(p0[3]));
    //return enlarge(b,Vec3fa(b.upper.w));
    //}

    /*! fill from list */
    __forceinline void fill(atomic_set<PrimRefBlockT<BezierPrim> >::block_iterator_unsafe& iter, Scene* scene, const bool list)
    {
      const BezierPrim& curve = *iter; iter++;
      const unsigned geomID = curve.geomID<0>();
      const unsigned primID = curve.primID<0>();
      const BezierCurves* in = (BezierCurves*) scene->get(geomID);
      //const Vec3fa& p0 = in->vertex(in->curve(primID),0);
      //const Vec3fa& p1 = in->vertex(in->curve(primID),1);
      const size_t vertexID = in->curve(primID);
      new (this) Bezier1iMB(vertexID,geomID,primID,list && !iter);
    }

  public:
    //const Vec3fa* p0;      //!< pointer to first control point (x,y,z,r) for time t0
    //const Vec3fa* p1;      //!< pointer to first control point (x,y,z,r) for time t1
    unsigned int vertexID; //!< index of start vertex
    unsigned int geom;  //!< geometry ID
    unsigned int prim;  //!< primitive ID
  };

  struct Bezier1iType : public PrimitiveType {
    Bezier1iType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };

  struct SceneBezier1i : public Bezier1iType
  {
    static SceneBezier1i type;
    BBox3fa update(char* prim, size_t num, void* geom) const;
  };

  struct Bezier1iMBType : public PrimitiveType {
    static Bezier1iMBType type;
    Bezier1iMBType ();
    size_t blocks(size_t x) const;
    size_t size(const char* This) const;
  };
}
