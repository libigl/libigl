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
#include "matte_textured_ispc.h"

namespace embree
{
  struct MatteTextured
  {
    static void* create(const Parms& parms)
    {
      ISPCRef Kd = parms.getTexture("Kd");
      const Vec2f s0 = parms.getVec2f("s0",Vec2f(0.0f,0.0f));
      const Vec2f ds = parms.getVec2f("ds",Vec2f(1.0f,1.0f));
      return ispc::MatteTextured__new(Kd.ptr,(ispc::vec2f&)s0,(ispc::vec2f&)ds);
    }
  };
}

