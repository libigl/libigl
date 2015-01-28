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

#include "scene_user_geometry.h"
#include "scene.h"

namespace embree
{
#if !defined (__MIC__)
  extern AccelSet::IntersectFunc4 ispcWrapperIntersect4;
  extern AccelSet::OccludedFunc4 ispcWrapperOccluded4;
#endif

#if defined(__TARGET_AVX__)
  extern AccelSet::IntersectFunc8 ispcWrapperIntersect8;
  extern AccelSet::OccludedFunc8 ispcWrapperOccluded8;
#endif

#if defined(__MIC__)
  extern AccelSet::IntersectFunc16 ispcWrapperIntersect16;
  extern AccelSet::OccludedFunc16 ispcWrapperOccluded16;
#endif

  UserGeometryScene::Base::Base (Scene* parent, GeometryTy ty, size_t items)
    : Geometry(parent,ty,1,RTC_GEOMETRY_STATIC), AccelSet(items)
  {
    enabling();
  }

  void UserGeometryScene::Base::enabling () { 
    atomic_add(&parent->numUserGeometries,+1); 
  }
  
  void UserGeometryScene::Base::disabling() { 
    atomic_add(&parent->numUserGeometries,-1); 
  }

  UserGeometryScene::UserGeometry::UserGeometry (Scene* parent, size_t items) 
    : Base(parent,USER_GEOMETRY,items),
      ispcPtr(NULL),
      ispcIntersect1(NULL), ispcIntersect4(NULL), ispcIntersect8(NULL), ispcIntersect16(NULL),
      ispcOccluded1(NULL), ispcOccluded4(NULL), ispcOccluded8(NULL), ispcOccluded16(NULL) {}
  
  void UserGeometryScene::UserGeometry::setUserData (void* ptr, bool ispc) 
  {
    if (ispc) ispcPtr = ptr;
    else intersectors.ptr = ptr;
    intersectors.boundsPtr = ptr;
  }
  
  void UserGeometryScene::UserGeometry::setBoundsFunction (RTCBoundsFunc bounds) {
    this->boundsFunc = bounds;
  }

  void UserGeometryScene::UserGeometry::setIntersectFunction (RTCIntersectFunc intersect1, bool ispc) {
    intersectors.intersector1.intersect = intersect1;
  }

  void UserGeometryScene::UserGeometry::setIntersectFunction4 (RTCIntersectFunc4 intersect4, bool ispc) 
  {
    if (ispc) {
#if !defined (__MIC__)
      intersectors.ptr = this;
      intersectors.intersector4.intersect = ispcWrapperIntersect4;
      ispcIntersect4 = (void*) intersect4;
#endif
    } else {
      intersectors.intersector4.intersect = intersect4;
    }
  }

  void UserGeometryScene::UserGeometry::setIntersectFunction8 (RTCIntersectFunc8 intersect8, bool ispc) 
  {
    if (ispc) {
#if defined(__TARGET_AVX__)
      intersectors.ptr = this;
      intersectors.intersector8.intersect = ispcWrapperIntersect8;
      ispcIntersect8 = (void*) intersect8;
#endif
    } else {
      intersectors.intersector8.intersect = intersect8;
    }
  }

  void UserGeometryScene::UserGeometry::setIntersectFunction16 (RTCIntersectFunc16 intersect16, bool ispc) 
  {
    if (ispc) {
#if defined(__MIC__)
      intersectors.ptr = this;
      intersectors.intersector16.intersect = ispcWrapperIntersect16;
      ispcIntersect16 = (void*) intersect16;
#endif
    } else {
      intersectors.intersector16.intersect = intersect16;
    }
  }

  void UserGeometryScene::UserGeometry::setOccludedFunction (RTCOccludedFunc occluded1, bool ispc) {
    intersectors.intersector1.occluded = occluded1;
  }

  void UserGeometryScene::UserGeometry::setOccludedFunction4 (RTCOccludedFunc4 occluded4, bool ispc) 
  {
    if (ispc) {
#if !defined (__MIC__)
      intersectors.ptr = this;
      intersectors.intersector4.occluded = ispcWrapperOccluded4;
      ispcOccluded4 = (void*) occluded4;
#endif
    } else {
      intersectors.intersector4.occluded = occluded4;
    }
  }

  void UserGeometryScene::UserGeometry::setOccludedFunction8 (RTCOccludedFunc8 occluded8, bool ispc) 
  {
    if (ispc) {
#if defined(__TARGET_AVX__)
      intersectors.ptr = this;
      intersectors.intersector8.occluded = ispcWrapperOccluded8;
      ispcOccluded8 = (void*) occluded8;
#endif
    } else {
      intersectors.intersector8.occluded = occluded8;
    }
  }

  void UserGeometryScene::UserGeometry::setOccludedFunction16 (RTCOccludedFunc16 occluded16, bool ispc) 
  {
    if (ispc) {
#if defined(__MIC__)
      intersectors.ptr = this;
      intersectors.intersector16.occluded = ispcWrapperOccluded16;
      ispcOccluded16 = (void*) occluded16;
#endif
    } else {
      intersectors.intersector16.occluded = occluded16;
    }
  }

  extern RTCBoundsFunc InstanceBoundsFunc;
  extern AccelSet::Intersector1 InstanceIntersector1;
  extern AccelSet::Intersector4 InstanceIntersector4;
  extern AccelSet::Intersector8 InstanceIntersector8;
  extern AccelSet::Intersector16 InstanceIntersector16;

  UserGeometryScene::Instance::Instance (Scene* parent, Accel* object) 
    : Base(parent,INSTANCES,1), local2world(one), world2local(one), object(object)
  {
    intersectors.ptr = this;
    intersectors.boundsPtr = this;
    boundsFunc = InstanceBoundsFunc;
    intersectors.intersector1 = InstanceIntersector1;
    intersectors.intersector4 = InstanceIntersector4; 
    intersectors.intersector8 = InstanceIntersector8; 
    intersectors.intersector16 = InstanceIntersector16;
  }
  
  void UserGeometryScene::Instance::setTransform(AffineSpace3f& xfm)
  {
    local2world = xfm;
    world2local = rcp(xfm);
  }
}
