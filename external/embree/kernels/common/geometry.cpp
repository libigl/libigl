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

#include "geometry.h"
#include "scene.h"

namespace embree
{
  Geometry::Geometry (Scene* parent, GeometryTy type, size_t numPrimitives, RTCGeometryFlags flags) 
    : parent(parent), type(type), numPrimitives(numPrimitives), id(0), flags(flags), state(ENABLING),
      intersectionFilter1(NULL), occlusionFilter1(NULL),
      intersectionFilter4(NULL), occlusionFilter4(NULL), ispcIntersectionFilter4(NULL), ispcOcclusionFilter4(NULL), 
      intersectionFilter8(NULL), occlusionFilter8(NULL), ispcIntersectionFilter8(NULL), ispcOcclusionFilter8(NULL), 
      intersectionFilter16(NULL), occlusionFilter16(NULL), ispcIntersectionFilter16(NULL), ispcOcclusionFilter16(NULL), 
      userPtr(NULL)
  {
    id = parent->add(this);
  }

  void Geometry::enable () 
  {
    if (parent->isStatic()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get enabled");
      return;
    }

    if (isDisabled()) {
      atomic_add(&parent->numIntersectionFilters4,(intersectionFilter4 != NULL) + (occlusionFilter4 != NULL));
      atomic_add(&parent->numIntersectionFilters8,(intersectionFilter8 != NULL) + (occlusionFilter8 != NULL));
      atomic_add(&parent->numIntersectionFilters16,(intersectionFilter16 != NULL) + (occlusionFilter16 != NULL));
    }

    switch (state) {
    case ENABLING:
      break;
    case ENABLED:
      break;
    case MODIFIED:
      break;
    case DISABLING: 
      state = MODIFIED;
      enabling();
      break;
    case DISABLED: 
      state = ENABLING;
      enabling();
      break;
    case ERASING:
      break;
    }
  }

  void Geometry::update() 
  {
    if (parent->isStatic()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get updated");
      return;
    }

    switch (state) {
    case ENABLING:
      break;
    case ENABLED:
      state = MODIFIED;
      break;
    case MODIFIED:
      break;
    case DISABLING: 
      break;
    case DISABLED: 
      break;
    case ERASING:
      break;
    }
  }

  void Geometry::disable () 
  {
    if (parent->isStatic()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get disabled");
      return;
    }

    if (isEnabled()) {
      atomic_sub(&parent->numIntersectionFilters4,(intersectionFilter4 != NULL) + (occlusionFilter4 != NULL));
      atomic_sub(&parent->numIntersectionFilters8,(intersectionFilter8 != NULL) + (occlusionFilter8 != NULL));
      atomic_sub(&parent->numIntersectionFilters16,(intersectionFilter16 != NULL) + (occlusionFilter16 != NULL));
    }

    switch (state) {
    case ENABLING:
      state = DISABLED;
      disabling();
      break;
    case ENABLED:
      state = DISABLING;
      disabling();
      break;
    case MODIFIED:
      state = DISABLING;
      disabling();
      break;
    case DISABLING: 
      break;
    case DISABLED: 
      break;
    case ERASING:
      break;
    }
  }

  void Geometry::erase () 
  {
    if (parent->isStatic()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get deleted");
      return;
    }

    switch (state) {
    case ENABLING:
      state = ERASING;
      disabling();
      break;
    case ENABLED:
      state = ERASING;
      disabling();
      break;
    case MODIFIED:
      state = ERASING;
      disabling();
      break;
    case DISABLING: 
      state = ERASING;
      break;
    case DISABLED: 
      state = ERASING;
      break;
    case ERASING:
      break;
    }
  }

  void Geometry::setIntersectionFilterFunction (RTCFilterFunc filter, bool ispc) 
  {
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    intersectionFilter1 = filter;
  }
    
  void Geometry::setIntersectionFilterFunction4 (RTCFilterFunc4 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters4,intersectionFilter4 != NULL);
    atomic_add(&parent->numIntersectionFilters4,filter != NULL);
    intersectionFilter4 = filter;
    if (ispc) ispcIntersectionFilter4 = (void*) filter; 
    else      ispcIntersectionFilter4 = NULL;
  }
    
  void Geometry::setIntersectionFilterFunction8 (RTCFilterFunc8 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters8,intersectionFilter8 != NULL);
    atomic_add(&parent->numIntersectionFilters8,filter != NULL);
    intersectionFilter8 = filter;
    if (ispc) ispcIntersectionFilter8 = (void*) filter; 
    else      ispcIntersectionFilter8 = NULL;
  }
  
  void Geometry::setIntersectionFilterFunction16 (RTCFilterFunc16 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters16,intersectionFilter16 != NULL);
    atomic_add(&parent->numIntersectionFilters16,filter != NULL);
    intersectionFilter16 = filter;
    if (ispc) ispcIntersectionFilter16 = (void*) filter; 
    else      ispcIntersectionFilter16 = NULL;
  }

  void Geometry::setOcclusionFilterFunction (RTCFilterFunc filter, bool ispc) 
  {
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    occlusionFilter1 = filter;
  }
    
  void Geometry::setOcclusionFilterFunction4 (RTCFilterFunc4 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters4,occlusionFilter4 != NULL);
    atomic_add(&parent->numIntersectionFilters4,filter != NULL);
    occlusionFilter4 = filter;
    if (ispc) ispcOcclusionFilter4 = (void*) filter; 
    else      ispcOcclusionFilter4 = NULL;
  }
    
  void Geometry::setOcclusionFilterFunction8 (RTCFilterFunc8 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters8,occlusionFilter8 != NULL);
    atomic_add(&parent->numIntersectionFilters8,filter != NULL);
    occlusionFilter8 = filter;
    if (ispc) ispcOcclusionFilter8 = (void*) filter; 
    else      ispcOcclusionFilter8 = NULL;
  }
  
  void Geometry::setOcclusionFilterFunction16 (RTCFilterFunc16 filter, bool ispc) 
  { 
    if (type != TRIANGLE_MESH && type != BEZIER_CURVES) {
      process_error(RTC_INVALID_OPERATION,"filter functions only supported for triangle meshes and hair geometries"); 
      return;
    }
    atomic_sub(&parent->numIntersectionFilters16,occlusionFilter16 != NULL);
    atomic_add(&parent->numIntersectionFilters16,filter != NULL);
    occlusionFilter16 = filter;
    if (ispc) ispcOcclusionFilter16 = (void*) filter; 
    else      ispcOcclusionFilter16 = NULL;
  }
}
