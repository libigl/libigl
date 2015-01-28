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

#include "ispc_wrapper_knc.h"

namespace embree
{
  typedef void (*ISPCIntersectFunc16)(void* ptr, RTCRay16& ray, size_t item, __mmask16 valid);
  typedef void (*ISPCOccludedFunc16 )(void* ptr, RTCRay16& ray, size_t item, __mmask16 valid);

  void ISPCWrapperKNC::intersect(const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay16& ray, size_t item) 
  {
    mic_i maski = *(mic_i*)valid;
    __mmask16 mask = maski != mic_i(0);
    ((ISPCIntersectFunc16)geom->ispcIntersect16)(geom->ispcPtr,ray,item,mask);
  }

  void ISPCWrapperKNC::occluded (const void* valid, const UserGeometryScene::UserGeometry* geom, RTCRay16& ray, size_t item) 
  {
    mic_i maski = *(mic_i*)valid;
    __mmask16 mask = maski != mic_i(0);
    ((ISPCOccludedFunc16)geom->ispcOccluded16)(geom->ispcPtr,ray,item,mask);
  }

  AccelSet::IntersectFunc16 ispcWrapperIntersect16 = (AccelSet::IntersectFunc16) ISPCWrapperKNC::intersect;
  AccelSet::OccludedFunc16  ispcWrapperOccluded16 =  (AccelSet::OccludedFunc16 ) ISPCWrapperKNC::occluded;
}
