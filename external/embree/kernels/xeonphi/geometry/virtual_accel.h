// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "primitive.h"

namespace embree
{
  struct AccelSetItem // FIXME: rename to Object
  {
  public:

    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& i, size_t end, Scene* scene) // FIXME: use nontemporal stores
    {
      const PrimRef& prim = prims[i];

      accel = (AccelSet*) (UserGeometryBase*) scene->get(prim.geomID());
      item  = prim.primID();

      i++;
    }

  public:
    AccelSet* accel;
    size_t item;
  };
}
