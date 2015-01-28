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

#ifndef __EMBREE_BUILD_SOURCE_H__
#define __EMBREE_BUILD_SOURCE_H__

#include "common/default.h"
#include "primref.h"

namespace embree
{
  /*! Interface of builder to geometry. */
  struct BuildSource
  {
    /*! return true if source is empty */
    virtual bool isEmpty () const = 0;

    /*! returns the number of groups */
    virtual size_t groups () const = 0;

    /*! returns the number of primitives per group */
    virtual size_t prims (size_t group, size_t* pNumVertices = NULL) const = 0;

    /*! calculates the bounding box of specified primitive of specified group */
    virtual const BBox3f bounds(size_t group, size_t prim) const = 0;

    /*! calculates the bounding box of specified primitive of specified group */
    virtual void bounds(size_t group, size_t begin, size_t end, BBox3f* bounds_o) const {}

    /*! splits a clipped primitive into two clipped primitives */
    virtual void split (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) const { 
      throw std::runtime_error("split not implemented"); 
    }

    /*! calculates number of primitives */
    size_t size() const 
    {
      size_t N = 0;
      size_t G = groups();
      for (size_t g=0; g<G; g++)
        N += prims(g);
      return N;
    }
  };

  __forceinline void splitTriangle(const PrimRef& prim, int dim, float pos, const Vec3fa& a, const Vec3fa& b, const Vec3fa& c, PrimRef& left_o, PrimRef& right_o)
  {
    BBox3f left = empty, right = empty;
    const Vec3fa v[3] = { a,b,c };

    /* clip triangle to left and right box by processing all edges */
    Vec3fa v1 = v[2];
    for (size_t i=0; i<3; i++)
    {
      Vec3fa v0 = v1; v1 = v[i];
      float v0d = v0[dim], v1d = v1[dim];
      
      if (v0d <= pos) left. extend(v0); // this point is on left side
      if (v0d >= pos) right.extend(v0); // this point is on right side

      if ((v0d < pos && pos < v1d) || (v1d < pos && pos < v0d)) // the edge crosses the splitting location
      {
        assert((v1d-v0d) != 0.0f);
        Vec3fa c = v0 + (pos-v0d)/(v1d-v0d)*(v1-v0);
        left.extend(c);
        right.extend(c);
      }
    }
    assert(!left.empty());  // happens if split does not hit triangle
    assert(!right.empty()); // happens if split does not hit triangle

    /* safe clip against current bounds */
    BBox3f bounds = prim.bounds();
    BBox3f cleft(min(max(left.lower,bounds.lower),bounds.upper),
                 max(min(left.upper,bounds.upper),bounds.lower));
    BBox3f cright(min(max(right.lower,bounds.lower),bounds.upper),
                  max(min(right.upper,bounds.upper),bounds.lower));

    new (&left_o ) PrimRef(cleft, prim.geomID(), prim.primID());
    new (&right_o) PrimRef(cright,prim.geomID(), prim.primID());
  }
}
#endif
