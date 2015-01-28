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

#include "ispc_wrapper_sse.h"

namespace embree
{
  typedef void (*ISPCIntersectFunc4)(void* ptr, RTCRay4& ray, size_t item, __m128 valid);
  typedef void (*ISPCOccludedFunc4 )(void* ptr, RTCRay4& ray, size_t item, __m128 valid);

  void ISPCWrapperSSE::intersect(const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay4& ray, size_t item) {
    assert(geom->ispcIntersect4);
    ((ISPCIntersectFunc4)geom->ispcIntersect4)(geom->ispcPtr,ray,item,*(__m128*)valid);
  }

  void ISPCWrapperSSE::occluded (const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay4& ray, size_t item) {
    assert(geom->ispcOccluded4);
    ((ISPCOccludedFunc4)geom->ispcOccluded4)(geom->ispcPtr,ray,item,*(__m128*)valid);
  }

  AccelSet::IntersectFunc4 ispcWrapperIntersect4 = (AccelSet::IntersectFunc4) ISPCWrapperSSE::intersect;
  AccelSet::OccludedFunc4  ispcWrapperOccluded4  = (AccelSet::OccludedFunc4 ) ISPCWrapperSSE::occluded;
}

