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

#include "scene.h"

namespace embree
{
  UserGeometryBase::UserGeometryBase (Scene* parent, GeometryTy ty, size_t items)
    : Geometry(parent,ty,1,RTC_GEOMETRY_STATIC), AccelSet(items)
  {
    enabling();
  }

  void UserGeometryBase::enabling () { 
    atomic_add(&parent->numUserGeometries1,numItems); 
  }
  
  void UserGeometryBase::disabling() { 
    atomic_add(&parent->numUserGeometries1,-(ssize_t)numItems); 
  }

  UserGeometry::UserGeometry (Scene* parent, size_t items) 
    : UserGeometryBase(parent,USER_GEOMETRY,items) {}
  
  void UserGeometry::setUserData (void* ptr, bool ispc) {
    intersectors.ptr = ptr;
  }
  
  void UserGeometry::setBoundsFunction (RTCBoundsFunc bounds) {
    this->boundsFunc = bounds;
  }

  void UserGeometry::setIntersectFunction (RTCIntersectFunc intersect1, bool ispc) {
    intersectors.intersector1.intersect = intersect1;
  }

  void UserGeometry::setIntersectFunction4 (RTCIntersectFunc4 intersect4, bool ispc) 
  {
    intersectors.intersector4.intersect = (void*)intersect4;
    intersectors.intersector4.ispc = ispc;
  }

  void UserGeometry::setIntersectFunction8 (RTCIntersectFunc8 intersect8, bool ispc) 
  {
    intersectors.intersector8.intersect = (void*)intersect8;
    intersectors.intersector8.ispc = ispc;
  }

  void UserGeometry::setIntersectFunction16 (RTCIntersectFunc16 intersect16, bool ispc) 
  {
    intersectors.intersector16.intersect = (void*)intersect16;
    intersectors.intersector16.ispc = ispc;
  }

  void UserGeometry::setOccludedFunction (RTCOccludedFunc occluded1, bool ispc) {
    intersectors.intersector1.occluded = occluded1;
  }

  void UserGeometry::setOccludedFunction4 (RTCOccludedFunc4 occluded4, bool ispc) 
  {
    intersectors.intersector4.occluded = (void*)occluded4;
    intersectors.intersector4.ispc = ispc;
  }

  void UserGeometry::setOccludedFunction8 (RTCOccludedFunc8 occluded8, bool ispc) 
  {
    intersectors.intersector8.occluded = (void*)occluded8;
    intersectors.intersector8.ispc = ispc;
  }

  void UserGeometry::setOccludedFunction16 (RTCOccludedFunc16 occluded16, bool ispc) 
  {
    intersectors.intersector16.occluded = (void*)occluded16;
    intersectors.intersector16.ispc = ispc;
  }

  extern RTCBoundsFunc InstanceBoundsFunc;
  extern AccelSet::Intersector1 InstanceIntersector1;
  extern AccelSet::Intersector4 InstanceIntersector4;
  extern AccelSet::Intersector8 InstanceIntersector8;
  extern AccelSet::Intersector16 InstanceIntersector16;

  Instance::Instance (Scene* parent, Accel* object) 
    : UserGeometryBase(parent,USER_GEOMETRY,1), local2world(one), world2local(one), object(object)
  {
    intersectors.ptr = this;
    boundsFunc = InstanceBoundsFunc;
    intersectors.intersector1 = InstanceIntersector1;
    intersectors.intersector4 = InstanceIntersector4; 
    intersectors.intersector8 = InstanceIntersector8; 
    intersectors.intersector16 = InstanceIntersector16;
  }
  
  void Instance::setTransform(AffineSpace3fa& xfm)
  {
    local2world = xfm;
    world2local = rcp(xfm);
  }
}
