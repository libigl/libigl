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
#include "pinholecamera_ispc.h"

namespace embree
{
  struct PinHoleCamera
  {
    static void* create(const Parms& parms)
    {
      const AffineSpace3f local2world = parms.getTransform("local2world");
      const float fov                 = parms.getFloat("angle",64.0f);
      const float aspectRatio         = parms.getFloat("aspectRatio",1.0f);
      return ispc::PinHoleCamera__new((ispc::vec3f&)local2world.l.vx,
                                      (ispc::vec3f&)local2world.l.vy,
                                      (ispc::vec3f&)local2world.l.vz,
                                      (ispc::vec3f&)local2world.p,
                                      fov,aspectRatio);
    }
  };
}

