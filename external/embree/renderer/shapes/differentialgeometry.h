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

#ifndef __EMBREE_DIFFERENTIAL_GEOMETRY_H__
#define __EMBREE_DIFFERENTIAL_GEOMETRY_H__

#include "rtcore/common/hit.h"

namespace embree
{
  /*! Contains additional shading information for hit points. */
  struct DifferentialGeometry : public Hit
  {
    /*! Default construction. */
    __forceinline DifferentialGeometry()
      : material(NULL), light(NULL) {}

    /*! Default construction. */
    __forceinline DifferentialGeometry(const Hit& hit)
      : Hit(hit), material(NULL), light(NULL) {}

  public:
    class Material*  material; //!< pointer to material of hit shape instance
    class AreaLight* light;    //!< pointer to area light of hit shape instance

  public:
    Vec3f P;         //!< Hit location in world coordinates.
    Vec3f Ng;        //!< Normalized geometric normal.
    Vec3f Ns;        //!< Normalized shading normal.
    Vec2f st;        //!< Hit location in surface parameter space.
    float error;     //!< Intersection error factor.
  };
}

#endif
