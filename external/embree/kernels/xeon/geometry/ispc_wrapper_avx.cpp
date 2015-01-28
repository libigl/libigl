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

#include "ispc_wrapper_avx.h"

namespace embree
{
  typedef void (*ISPCIntersectFunc8)(void* ptr, RTCRay8& ray, size_t item, __m256 valid);
  typedef void (*ISPCOccludedFunc8 )(void* ptr, RTCRay8& ray, size_t item, __m256 valid);

  void ISPCWrapperAVX::intersect(const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay8& ray, size_t item) {
    assert(geom->ispcIntersect8);
    ((ISPCIntersectFunc8)geom->ispcIntersect8)(geom->ispcPtr,ray,item,*(__m256*)valid);
  }

  void ISPCWrapperAVX::occluded (const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay8& ray, size_t item) {
    assert(geom->ispcOccluded8);
    ((ISPCOccludedFunc8)geom->ispcOccluded8)(geom->ispcPtr,ray,item,*(__m256*)valid);
  }

  AccelSet::IntersectFunc8 ispcWrapperIntersect8 = (AccelSet::IntersectFunc8) ISPCWrapperAVX::intersect;
  AccelSet::OccludedFunc8  ispcWrapperOccluded8  = (AccelSet::OccludedFunc8 ) ISPCWrapperAVX::occluded;
}
