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
#include "obj_ispc.h"

namespace embree
{
  struct Obj
  {
    static void* create(const Parms& parms)
    {
      /*! opacity texture and coefficient */
      ISPCRef map_d = parms.getTexture("map_d");  
      const float d = parms.getFloat("d", 1.0f);

      /*! diffuse reflectance texture and coefficient */
      ISPCRef map_Kd = parms.getTexture("map_Kd");  
      const Color Kd = parms.getColor("Kd", one);

      /*! specular reflectance texture and coefficient */
      ISPCRef map_Ks = parms.getTexture("map_Ks");  
      const Color Ks = parms.getColor("Ks", zero);

      /*! specular exponent texture and coefficient */
      ISPCRef map_Ns = parms.getTexture("map_Ns");  
      const float Ns = parms.getFloat("Ns", 10.0f);

      /*! get bump map */
      ISPCRef map_Bump = parms.getTexture("map_Bump");

      return ispc::Obj__new(map_d.ptr,d,
                            map_Kd.ptr,(ispc::vec3f&)Kd,
                            map_Ks.ptr,(ispc::vec3f&)Ks,
                            map_Ns.ptr,Ns,
                            map_Bump.ptr);
    }
  };
}

