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

#ifndef __EMBREE_INTERSECTOR8_H__
#define __EMBREE_INTERSECTOR8_H__

#include <stddef.h>
#include <immintrin.h>

namespace embree
{
  class Accel;
  struct Ray8;

  /*! Packet of 8 rays interface to the traverser. A closest intersection
   *  point of a ray with the geometry can be found. A ray can also be
   *  tested for occlusion by any geometry. */
  class Intersector8 {
  public:

    /*! name for this interface */
    static const char* const name;

    /*! Intersect function pointer type. */
    typedef void (*intersectFunc)(const Intersector8* This, /*!< pointer to this intersector */ 
                                  Ray8& ray,          /*!< Ray to shoot. */
                                  const __m256 valid        /*!< determines active rays */);
    
    /*! Occluded function pointer type. */
    typedef __m256 (*occludedFunc) (const Intersector8* This, /*!< pointer to this intersector */ 
                                    Ray8& ray,          /*!< Ray to test occlusion for. */
                                    const __m256 valid        /*!< determines active rays */);

  public:
  
    /*! Default constructor. */
    __forceinline Intersector8 ()
      : intersectPtr(NULL), occludedPtr(NULL) {}

    /*! Constructs the intersector from two function pointers. */
    __forceinline Intersector8 (intersectFunc intersect, occludedFunc occluded)
      : intersectPtr(intersect), occludedPtr(occluded) {}

    /*! Intersects the ray packet with the geometry and returns the hit
     *  information. */
    __forceinline void intersect(const __m256& valid, Ray8& ray) const {
      intersectPtr(this,ray,valid);
    }

    /*! Tests the ray packet for occlusion with the scene. */
    __forceinline __m256 occluded (const __m256& valid, Ray8& ray) const {
      return occludedPtr(this,ray,valid);
    }

  public:
    intersectFunc intersectPtr;  /*!< Pointer to intersect function */
    occludedFunc  occludedPtr;   /*!< Pointer to occluded function */
  };

  typedef Intersector8 RTCIntersector8;
}

#endif
