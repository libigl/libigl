// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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
  template class Base<HeuristicBinning<0> >;            \
  template class Base<HeuristicBinning<2> >;         \
  template class Base<HeuristicBinning<3> >;         \
  template class Base<HeuristicSpatial<0> >;         \
  template class Base<HeuristicSpatial<2> >;         \
  template class Base<HeuristicSpatial<3> >;

#define INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(Splitter)     \
  template class Splitter<HeuristicBinning<0>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicBinning<2>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicBinning<3>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<0>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<2>, atomic_set<PrimRefBlock> >; \
  template class Splitter<HeuristicSpatial<3>, atomic_set<PrimRefBlock> >;

#define INSTANTIATE_TEMPLATE_BUILDER_HELP(Builder,Heuristic)            \
  template class Builder<Heuristic<Triangle1i  ::logBlockSize>,Triangle1i>; \
  template class Builder<Heuristic<Triangle4i  ::logBlockSize>,Triangle4i>; \
  template class Builder<Heuristic<Triangle1   ::logBlockSize>,Triangle1>; \
  template class Builder<Heuristic<Triangle4   ::logBlockSize>,Triangle4>; \
  template class Builder<Heuristic<Triangle8   ::logBlockSize>,Triangle8>; \
  template class Builder<Heuristic<Triangle1Xfm::logBlockSize>,Triangle1Xfm>; \
  template class Builder<Heuristic<Triangle4Xfm::logBlockSize>,Triangle4Xfm>;

#define INSTANTIATE_TEMPLATE_BUILDER(Builder)                   \
  INSTANTIATE_TEMPLATE_BUILDER_HELP(Builder,HeuristicBinning)   \
  INSTANTIATE_TEMPLATE_BUILDER_HELP(Builder,HeuristicSpatial)

#define INSTANTIATE_TEMPLATE_BUILDER_HELP_MB(Builder,Heuristic)         \
  template class Builder<Heuristic<Triangle4i::logBlockSize>,Triangle4i>;

#define INSTANTIATE_TEMPLATE_BUILDER_MB(Builder)                        \
  INSTANTIATE_TEMPLATE_BUILDER_HELP_MB(Builder,HeuristicBinning)

#endif
