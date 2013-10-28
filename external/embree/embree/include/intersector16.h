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

#ifndef __EMBREE_INTERSECTOR16_H__
#define __EMBREE_INTERSECTOR16_H__

#include <stddef.h>
#include <immintrin.h>
#include <zmmintrin.h>

namespace embree
{
  class Accel;
  struct Ray16;

  /*! Packet of 16 rays interface to the traverser. A closest intersection
   *  point of a ray with the geometry can be found. A ray can also be
   *  tested for occlusion by any geometry. */
  class Intersector16 {
  public:
    
    /*! name for this interface */
    static const char* const name;

    /*! Intersects the ray packet with the geometry and returns the hit
     *  information. */
    typedef void (*intersectFunc)(const Intersector16* This, /*!< pointer to this intersector */
                                  Ray16& ray,                /*!< Ray to shoot. */
                                  const __mmask valid       /*!< determines active rays */);
    
    /*! Tests the ray packet for occlusion with the scene. \returns occluded mask */
    typedef __mmask (*occludedFunc) (const Intersector16* This, /*!< pointer to this intersector */ 
                                     Ray16& ray,                /*!< Ray to test occlusion for. */
                                     const __mmask valid       /*!< determines active rays */);

  public:
  
    /*! Default constructor. */
    __forceinline Intersector16 ()
      : intersectPtr(NULL), occludedPtr(NULL)  {}

    /*! Constructs the intersector from two function pointers. */
    __forceinline Intersector16 (intersectFunc intersect, occludedFunc occluded)
      : intersectPtr(intersect), occludedPtr(occluded) {}

    /*! Intersects the ray packet with the geometry and returns the hit
     *  information. */
    __forceinline void intersect(const __mmask& valid, Ray16& ray) const {
      intersectPtr(this,ray,valid);
    }

    /*! Tests the ray packet for occlusion with the scene. */
    __forceinline __mmask occluded (const __mmask& valid, Ray16& ray) const {
      return occludedPtr(this,ray,valid);
    }

  public:
    intersectFunc intersectPtr;  /*!< Pointer to intersect function */
    occludedFunc  occludedPtr;   /*!< Pointer to occluded function */
  };

  typedef Intersector16 RTCIntersector16;
}

#endif
