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
#include "metallicpaint_ispc.h"

namespace embree
{
  struct MetallicPaint
  {
    static void* create(const Parms& parms)
    {
      const Color shadeColor    = parms.getColor("shadeColor",one);
      const Color glitterColor  = parms.getColor("glitterColor",zero);
      const float glitterSpread = parms.getFloat("glitterSpread",1.0f);
      const float eta           = parms.getFloat("eta",1.4f);
      return ispc::MetallicPaint__new((ispc::vec3f&)shadeColor,
                                      (ispc::vec3f&)glitterColor,
                                      glitterSpread,
                                      eta);
    }
  };
}

