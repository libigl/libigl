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

#ifndef __EMBREE_DIFFERENTIAL_GEOMETRY_H__
#define __EMBREE_DIFFERENTIAL_GEOMETRY_H__

#include "default.h" 

namespace embree
{
  typedef int light_mask_t;

  /*! Contains additional shading information for hit points. */
  struct DifferentialGeometry
  {
    /*! Default construction. */
    __forceinline DifferentialGeometry()
      : material(NULL), light(NULL) {}

  public:
    class Material*  material; //!< pointer to material of hit shape instance
    class AreaLight* light;    //!< pointer to area light of hit shape instance

  public:
    Vector3f P;                //!< Hit location in world coordinates.
    Vector3f Tx;               //!< Tangent in x direction
    Vector3f Ty;               //!< Tangent in y direction
    Vector3f Ng;               //!< Normalized geometry normal.
    mutable Vector3f Ns;       //!< Normalized shading normal.
    Vec2f st;                  //!< Hit location in surface parameter space.
    float error;               //!< Intersection error factor.
    light_mask_t illumMask;    //!< bit mask which light we're interested in
    light_mask_t shadowMask;   //!< bit mask which light we're interested in
  };
}

#endif
