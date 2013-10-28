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

#ifndef __EMBREE_INTERSECTOR1_H__
#define __EMBREE_INTERSECTOR1_H__

#include <stddef.h>

namespace embree
{
  class Accel;
  struct Ray;

  /*! Single ray interface to the traverser. A closest intersection
   *  point of a ray with the geometry can be found. A ray can also be
   *  tested for occlusion by any geometry. */
  class Intersector1 {
  public:

    /*! name for this interface */
    static const char* const name;

    /*! Intersect function pointer type. */
    typedef void (*intersectFunc)(const Intersector1* This,  /*!< pointer to this intersector */
                                  Ray& ray             /*!< Ray to shoot. */);
    
    /*! Occluded function pointer type. */
    typedef bool (*occludedFunc) (const Intersector1* This, /*!< pointer to this intersector */ 
                                  Ray& ray            /*!< Ray to test occlusion for. */);

  public:

    /*! Default constructor. */
    __forceinline Intersector1 ()
      : intersectPtr(NULL), occludedPtr(NULL) {}

    /*! Constructs the intersector from two function pointers. */
    __forceinline Intersector1 (intersectFunc intersect, occludedFunc occluded)
      : intersectPtr(intersect), occludedPtr(occluded) {}

    /*! Intersects the ray with the geometry and returns the hit
     *  information. */
    __forceinline void intersect(Ray& ray) const {
      intersectPtr(this,ray);
    }

    /*! Tests the ray for occlusion with the scene. */
    __forceinline bool occluded (Ray& ray) const {
      return occludedPtr(this,ray);
    }

  public:
    intersectFunc intersectPtr;  /*!< Pointer to intersect function */
    occludedFunc  occludedPtr;   /*!< Pointer to occluded function */
  };

  typedef Intersector1 RTCIntersector1;
}

#endif
