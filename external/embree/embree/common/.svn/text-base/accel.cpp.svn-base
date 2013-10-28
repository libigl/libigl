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

#include "accel.h"
#include "../geometry/geometry.h"

namespace embree
{
  Accel::Accel (RTCGeometry* geom, const TriangleType& trity)
    : geom(geom), trity(trity), vertices(NULL), numVertices(0)
  {
    if (trity.needVertices) {
      this->vertices = (Vec3fa*) geom->getVertices();
      this->numVertices = geom->getNumVertices();
    }
  }
}
