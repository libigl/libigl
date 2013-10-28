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

#pragma once

#include "api/parms.h"
#include "velvet_ispc.h"

namespace embree
{
  struct Velvet
  {
    static void* create(const Parms& parms)
    {
      const Color reflectance = parms.getColor("reflectance",one);
      const Color horizonScatteringColor = parms.getColor("horizonScatteringColor",one);
      const float horizonScatteringFallOff = parms.getFloat("horizonScatteringFallOff",zero);
      const float backScattering = parms.getFloat("backScattering",zero);
      return ispc::Velvet__new((ispc::vec3f&)reflectance,
                               (ispc::vec3f&)horizonScatteringColor,
                               horizonScatteringFallOff,
                               backScattering);
    }
  };
}

