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

#ifndef __EMBREE_HEURISTICS_H__
#define __EMBREE_HEURISTICS_H__

#include "primref.h"
#include "primrefblock.h"

#include "heuristic_binning.h"
#include "heuristic_spatial.h"

#define INSTANTIATE_TEMPLATE_BY_HEURISTIC(Base)         \
  template class Base<HeuristicBinning<0> >;                   \
  template class Base<HeuristicBinning<2> >;                   \
  template class Base<HeuristicBinning<3> >;                   \
  template class Base<HeuristicSpatial<0> >;                   \
  template class Base<HeuristicSpatial<2> >;                   \
  template class Base<HeuristicSpatial<3> >;

#define INSTANTIATE_TEMPLATE_BY_BVH_AND_HEURISTIC(Base,BVH)         \
  template class Base<BVH, HeuristicBinning<0> >;                   \
  template class Base<BVH, HeuristicBinning<2> >;                   \
  template class Base<BVH, HeuristicBinning<3> >;                   \
  template class Base<BVH, HeuristicSpatial<0> >;                   \
  template class Base<BVH, HeuristicSpatial<2> >;                   \
  template class Base<BVH, HeuristicSpatial<3> >;

#define INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(Splitter)     \
  template class Splitter<HeuristicBinning<0>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicBinning<2>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicBinning<3>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<0>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<2>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<3>, atomic_set<PrimRefBlock> >;

#endif
