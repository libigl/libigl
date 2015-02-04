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

#include "common/subdiv/subdivpatch1base.h"
#include "kernels/xeon/geometry/primitive.h"

namespace embree
{

  struct __aligned(64) SubdivPatch1Cached : public SubdivPatch1Base
  {

  SubdivPatch1Cached (const CatmullClarkPatch& ipatch,
                      const unsigned int gID,
                      const unsigned int pID,
                      const SubdivMesh *const mesh,
                      const Vec2f uv[4],
                      const float edge_level[4]) : SubdivPatch1Base(ipatch,gID,pID,mesh,uv,edge_level)
      {
      }

    struct Type : public PrimitiveType 
    {
      Type ();
      size_t blocks(size_t x) const; 
      size_t size(const char* This) const;
    };

    static Type type;

    /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return N; }

  };

};
