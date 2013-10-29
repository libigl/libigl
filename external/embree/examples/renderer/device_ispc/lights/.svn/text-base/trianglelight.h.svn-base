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
#include "trianglelight_ispc.h"

namespace embree
{
  struct TriangleLight
  {
    static void* create(const Parms& parms)
    {
      const Vector3f v0 = parms.getVector3f("v0");
      const Vector3f v1 = parms.getVector3f("v1");
      const Vector3f v2 = parms.getVector3f("v2");
      const Color L = parms.getColor("L");
      return ispc::TriangleLight__new((ispc::vec3f&)v0,(ispc::vec3f&)v1,(ispc::vec3f&)v2,(ispc::vec3f&)L);
    }
  };
}

