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

#include "subdivpatch1.h"
#include "common/scene.h"

namespace embree
{
  SubdivPatch1::Type SubdivPatch1::type;
  
  SubdivPatch1::Type::Type () 
    : PrimitiveType("subdivpatch1",sizeof(SubdivPatch1),1,false,1) {} 
  
  size_t SubdivPatch1::Type::blocks(size_t x) const {
    return x;
  }
    
  size_t SubdivPatch1::Type::size(const char* This) const {
    return 1;
  }
}
