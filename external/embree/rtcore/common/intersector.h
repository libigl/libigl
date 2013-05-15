// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_INTERSECTOR_H__
#define __EMBREE_INTERSECTOR_H__

#include "ray.h"
#include "hit.h"

namespace embree
{
  /*! Single ray interface to the traverser. A closest intersection
   *  point of a ray with the geometry can be found. A ray can also be
   *  tested for occlusion by any geometry. */
  class Intersector : public RefCount {
  public:

    /*! name for this interface */
    static const char* const name;

    /*! A virtual destructor is required. */
    virtual ~Intersector() {}

    /*! Intersects the ray with the geometry and returns the hit
     *  information. */
    virtual void intersect(const Ray& ray,   /*!< Ray to shoot. */
                           Hit& hit          /*!< Hit result.   */) const = 0;

    /*! Tests the ray for occlusion with the scene. */
    virtual bool occluded (const Ray& ray    /*!< Ray to test occlusion for. */) const = 0;
  };
}

#endif
