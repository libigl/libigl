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

#ifndef __EMBREE_ACCEL_H__
#define __EMBREE_ACCEL_H__

#include "registry_builder.h"
#include "registry_intersector.h"
#include "../geometry/triangle.h"

namespace embree
{
  struct RTCGeometry;

  /*! Abstract acceleration structure interface. */
  class Accel : public RefCount {
    ALIGNED_CLASS;
  public:
    
    /*! Default constructor. */
    Accel (RTCGeometry* geom, const TriangleType& trity);

  public:

    /*! return name of this acceleration structure */
    virtual const std::string name() const = 0;

    /*! prints statistics */
    virtual void print() = 0;
    
    __forceinline       void* nodePtr()       { return NULL; }
    __forceinline const void* nodePtr() const { return NULL; }

    __forceinline       void* triPtr()       { return NULL; }
    __forceinline const void* triPtr() const { return NULL; }

  public:
    RTCGeometry* geom;             //!< geometry this acceleration structure is build over
    const TriangleType& trity;     //!< triangle type stored in BVH
    const Vec3fa* vertices;        //!< Pointer to vertex array.
    size_t numVertices;            //!< Number of vertices
    BBox3f bounds;                 //!< Bounding box of geometry.
  };
}

#endif
