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

#ifndef __EMBREE_SPLITTER_H__
#define __EMBREE_SPLITTER_H__

#include "../geometry/geometry.h"
#include "primrefalloc.h"
#include "primrefblock.h"

namespace embree
{
  /*! Splits a list of build primitives into two lists using a single thread. */
  template<typename Heuristic>      
    class Splitter
  {
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
  public:
    /*! the constructor directly performs the split */
    Splitter (size_t thread, PrimRefAlloc* alloc, const RTCGeometry* geom,
              atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& psplit);

    /*! perform specified split */
    void split(size_t thread, PrimRefAlloc* alloc, const RTCGeometry* geom,
               atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

    /*! perform specified spatial split */
    void split_spatial(size_t thread, PrimRefAlloc* alloc, const RTCGeometry* geom,
                       atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);
    
    /* output data */
  public:
    atomic_set<PrimRefBlock> lprims;  //!< left primitive list
    atomic_set<PrimRefBlock> rprims;  //!< right primitive list
    PrimInfo linfo;                   //!< left bounding information of primitives
    PrimInfo rinfo;                   //!< right bounding information of primitives
    Split lsplit;                     //!< best split for left primitives
    Split rsplit;                     //!< best split for right primitives
  };
}

#endif
