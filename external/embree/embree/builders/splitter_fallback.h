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

#ifndef __EMBREE_SPLITTER_FALLBACK_H__
#define __EMBREE_SPLITTER_FALLBACK_H__

#include "../geometry/geometry.h"
#include "primrefalloc.h"
#include "primrefblock.h"

namespace embree
{
  /*! Splits a list of build primitives into two arbitrary lists. */
  template<typename Heuristic, typename PrimRefBlockList>      
    class FallBackSplitter
  {
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
  public:

    /*! enforce some object median split */
    static void split(size_t threadIndex, PrimRefAlloc* alloc, const RTCGeometry* geom,
                      PrimRefBlockList& prims, const PrimInfo& pinfo,
                      PrimRefBlockList& lprims_o, PrimInfo& linfo_o, Split& lsplit_o,
                      PrimRefBlockList& rprims_o, PrimInfo& rinfo_o, Split& rsplit_o);

    /*! enforce some object median split */
    static void split(size_t threadIndex, PrimRefAlloc* alloc, const RTCGeometry* geom,
                      PrimRefBlockList& prims, const PrimInfo& pinfo,
                      PrimRefBlockList& lprims_o, PrimInfo& linfo_o, 
                      PrimRefBlockList& rprims_o, PrimInfo& rinfo_o);
  };
}

#endif
