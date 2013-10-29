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

#ifndef __EMBREE_ACCEL_VIRTUAL_OBJECT_H__
#define __EMBREE_ACCEL_VIRTUAL_OBJECT_H__

#include "triangle.h"
#include "virtual_scene.h"
#include "../common/accel.h"

namespace embree
{
  /*! Object with transformation. */
  struct VirtualObject
  {
     /*! Name of intersector */
    static const char* name() { return "virtual"; }

    /*! block size */
    static const size_t blockSize = 1;
    static const size_t logBlockSize = 0;

    /*! Determines if we need the vertex array. */
    static const bool needVertices = true;

    /*! Cost of ray/triangle intersection. */
    static const int intCost = 1;

    /*! virtual interface to query information about the object type */
    static const struct Type : public TriangleType
    {
      Type () : TriangleType("virtual",sizeof(VirtualObject),1,true,1) {}

      size_t blocks(size_t x) const {
        return x;
      }
      
      size_t size(const char* This) const {
        return ((VirtualObject*)This)->size();
      }
      
      float area(const char* This, const RTCGeometry* geom) const {
        return ((VirtualObject*)This)->area(geom);
      }
      
      void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const {
        ((VirtualObject*)This)->pack(prims,(const VirtualScene*)geom);
      }
    } type;

  public:

    /*! Default constructor. */
    __forceinline VirtualObject () : id(0) {}

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const {
      return 1;
    }

    /*! Computes the area of the triangle. */
    __forceinline float area(const RTCGeometry* geom_i) {
      const VirtualScene* geom = (const VirtualScene*) geom_i;
      return embree::area(geom->get(id).worldBounds);
    }

    /*! Packs triangle taken from primitive list. */
    template<typename Iterator>
    __forceinline void pack(Iterator& prims, const VirtualScene* geom) {
      const PrimRef& prim = *prims; prims++; id = (int) prim.id();
    }

  public:
    int id; // index of object
  };
}

#endif
