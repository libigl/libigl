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

#ifndef __EMBREE_ACCEL_TRIANGLE_H__
#define __EMBREE_ACCEL_TRIANGLE_H__

#include "../common/default.h"
#include "../include/embree.h"

#include "../common/ray.h"

#if defined(__SSE__)
#include "../common/ray4.h"
#endif

#if defined(__AVX__)
#include "../common/ray8.h"
#endif

#if defined(__MIC__)
#include "../common/ray16.h"
#endif

#include "../builders/primrefblock.h"

namespace embree
{
  struct RTCGeometry;

  struct TriangleType
  {
    /*! constructs the triangle type */
    TriangleType (const char* name, size_t bytes, size_t blockSize, bool needVertices, int intCost) 
      : name(name), bytes(bytes), blockSize(blockSize), needVertices(needVertices), intCost(intCost) {}

    /*! Computes the number of blocks required to store a number of triangles. */
    virtual size_t blocks(size_t x) const = 0;

    /*! Returns the number of stored triangles. */
    virtual size_t size(const char* This) const = 0;

    /*! Computes the area of the triangle. */
    virtual float area(const char* This, const RTCGeometry* geom) const = 0;

    /*! Packs triangle taken from primitive list. */
    virtual void pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, const RTCGeometry* geom) const = 0;

    /*! Computes the bounds of timestep t0 and t1 of a moving triangle. */
    virtual std::pair<BBox3f,BBox3f> bounds(char* This, const RTCGeometry* geom) const { return std::pair<BBox3f,BBox3f>(empty,empty); }

    std::string name;       //!< name of this triangle type
    size_t bytes;           //!< number of bytes of the triangle data
    size_t blockSize;       //!< block size
    bool   needVertices;    //!< determines if we need the vertex array
    int    intCost;         //!< cost of one ray/triangle intersection
  };
}

#endif


