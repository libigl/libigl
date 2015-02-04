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

#include "geometry/bezier1v.h"
#include "builders/primrefalloc.h"
#include "priminfo.h"

namespace embree
{
  namespace isa
  {
    /*! Performs fallback splits */
    struct FallBackSplit
    {
      typedef atomic_set<PrimRefBlockT<PrimRef> > PrimRefList;
      typedef atomic_set<PrimRefBlockT<BezierPrim> > BezierRefList;
      
      __forceinline FallBackSplit (const NAABBox3fa& bounds0, size_t num0, const NAABBox3fa& bounds1, size_t num1)
	: bounds0(bounds0), num0(num0), bounds1(bounds1), num1(num1) {}
      
      /*! finds some partitioning */
      static FallBackSplit find(size_t threadIndex, PrimRefBlockAlloc<BezierPrim>& alloc, 
				BezierRefList& prims, BezierRefList& lprims_o, PrimInfo& linfo_o, BezierRefList& rprims_o, PrimInfo& rinfo_o);
      
      /*! finds some partitioning */
      static FallBackSplit find(size_t threadIndex, PrimRefBlockAlloc<PrimRef>& alloc, 
				PrimRefList& prims, PrimRefList& lprims_o, PrimInfo& linfo_o, PrimRefList& rprims_o, PrimInfo& rinfo_o);
      
    public:
      size_t num0, num1;
      NAABBox3fa bounds0, bounds1;
    };
  }
}
