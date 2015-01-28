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

#ifndef __EMBREE_BUILDER_H__
#define __EMBREE_BUILDER_H__

#include "common/default.h"

namespace embree
{
  class Builder : public RefCount {
  public:
    Builder () : needAllThreads(false) {}
    virtual void build(size_t threadIndex, size_t threadCount) = 0;
  public:
    bool needAllThreads;
  };

#define ADD_BUILDER(NAME,BUILDER,LEAFMIN,LEAFMAX)              \
  builders.add(ISA,NAME,BUILDER,LEAFMIN,LEAFMAX);

#define DECLARE_TRIANGLEMESH_BUILDER(symbol)                            \
  namespace isa   { extern Builder* symbol(void* accel, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace sse41 { extern Builder* symbol(void* accel, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace avx   { extern Builder* symbol(void* accel, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace avx2  { extern Builder* symbol(void* accel, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize); } \
  void symbol##_error() { std::cerr << "Error: builder " << TOSTRING(symbol) << " not supported no your CPU" << std::endl; } \
  TriangleMeshBuilderFunc symbol = (TriangleMeshBuilderFunc) symbol##_error;
  
#define DECLARE_BUILDER(symbol)                                         \
  namespace isa   { extern Builder* symbol(void* accel, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace sse41 { extern Builder* symbol(void* accel, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace avx   { extern Builder* symbol(void* accel, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize); } \
  namespace avx2  { extern Builder* symbol(void* accel, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize); } \
  void symbol##_error() { std::cerr << "Error: builder " << TOSTRING(symbol) << " not supported no your CPU" << std::endl; } \
  BuilderFunc symbol = (BuilderFunc) symbol##_error;
}

#endif

