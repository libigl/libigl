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
#include "common/subdiv/bspline_patch.h"
#include "common/subdiv/gregory_patch.h"

#define SUBDIVISION_LEVEL 3 // FIXME: removes

namespace embree
{
  struct SubdivPatch1
  {
    struct Type : public PrimitiveType 
    {
      Type ();
      size_t blocks(size_t x) const; 
      size_t size(const char* This) const;
    };
    
    static Type type;
    
  public:
    
    enum {
      REGULAR_PATCH = 1,
      HAS_BORDERS   = 2,
      HAS_CREASES   = 4
    };
    
    /*! Default constructor. */
    //__forceinline SubdivPatch1 () {}
    
    /*! Construction from vertices and IDs. */
    __forceinline SubdivPatch1 (const SubdivMesh::HalfEdge* edge, 
                                const BufferT<Vec3fa>& vertices, 
                                const unsigned int geom, 
                                const unsigned int prim, 
                                const unsigned int subdivision_level,
                                const bool last)
      : first_half_edge(edge), vertices(vertices), geom(geom), prim(prim | (last << 31)), subdivision_level(subdivision_level)
    {
      flags = 0;
      if (first_half_edge->isRegularFace()) 
      {
        flags |= REGULAR_PATCH;
#if 1
        BSplinePatch patch;
        
        init( patch );
#endif
      }     
    }
    
    __forceinline bool isRegular() const
    {
      return (flags & REGULAR_PATCH) == REGULAR_PATCH;
    }
    
    
    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return N; }
    
    /*! return geometry ID */
    template<bool list>
    __forceinline unsigned int geomID() const { 
      return geom; 
    }
    
    /*! return primitive ID */
    template<bool list>
    __forceinline unsigned int primID() const { 
      if (list) return prim & 0x7FFFFFFF; 
      else      return prim; 
    }
    
    /*! checks if this is the last primitive in list leaf mode */
    __forceinline int last() const { 
      return prim & 0x80000000; 
    }
    
    /*! builder interface to fill primitive */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      const PrimRef& prim = *prims;
      prims++;
      
      const unsigned int last   = list && !prims;
      const unsigned int geomID = prim.geomID();
      const unsigned int primID = prim.primID();
      const SubdivMesh* const subdiv_mesh = scene->getSubdivMesh(geomID);
      new (this) SubdivPatch1(subdiv_mesh->getHalfEdge(primID),
                              subdiv_mesh->getVertexBuffer(),
                              geomID,
                              primID,
                              SUBDIVISION_LEVEL,
                              last); 
    }
    
    /*! builder interface to fill primitive */
    __forceinline void fill(const PrimRef* prims, size_t& i, size_t end, Scene* scene, const bool list)
    {
      const PrimRef& prim = prims[i];
      i++;
      
      const unsigned int last = list && i >= end;
      const unsigned int geomID = prim.geomID();
      const unsigned int primID = prim.primID();
      const SubdivMesh* const subdiv_mesh = scene->getSubdivMesh(geomID);
      new (this) SubdivPatch1(subdiv_mesh->getHalfEdge(primID),
                              subdiv_mesh->getVertexBuffer(),
                              geomID,
                              primID,
                              SUBDIVISION_LEVEL,
                              last); 
    }
    
    __forceinline void init( CatmullClarkPatch& patch) const
    {
      for (size_t i=0; i<4; i++)
        patch.ring[i].init(first_half_edge + i,vertices);
    }
    
    __forceinline void init( BSplinePatch& cc_patch) const
    {
      cc_patch.init(first_half_edge,vertices);
    }
    
  public:
    const SubdivMesh::HalfEdge* first_half_edge;  //!< pointer to first half edge of this patch
    const BufferT<Vec3fa>& vertices;                       //!< pointer to vertex array
    unsigned int subdivision_level;
    unsigned int flags;
    unsigned int geom;                            //!< geometry ID of the subdivision mesh this patch belongs to
    unsigned int prim;                            //!< primitive ID of this subdivision patch
  };
}
