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

#ifndef __EMBREE_USER_GEOMETRY_SCENE_H__
#define __EMBREE_USER_GEOMETRY_SCENE_H__

#include "common/default.h"
#include "common/accel.h"
#include "common/accelset.h"
#include "common/geometry.h"

namespace embree
{
  namespace UserGeometryScene
  {
    struct Base : public Geometry, public AccelSet
    {

    public:
      Base (Scene* parent, GeometryTy ty, size_t items);
      
      void enabling ();
      void disabling();
    };

    struct UserGeometry : public Base
    {
    public:
      UserGeometry (Scene* parent, size_t items); 
      virtual void setUserData (void* ptr, bool ispc);
      virtual void setBoundsFunction (RTCBoundsFunc bounds);
      virtual void setIntersectFunction (RTCIntersectFunc intersect, bool ispc);
      virtual void setIntersectFunction4 (RTCIntersectFunc4 intersect4, bool ispc);
      virtual void setIntersectFunction8 (RTCIntersectFunc8 intersect8, bool ispc);
      virtual void setIntersectFunction16 (RTCIntersectFunc16 intersect16, bool ispc);
      virtual void setOccludedFunction (RTCOccludedFunc occluded, bool ispc);
      virtual void setOccludedFunction4 (RTCOccludedFunc4 occluded4, bool ispc);
      virtual void setOccludedFunction8 (RTCOccludedFunc8 occluded8, bool ispc);
      virtual void setOccludedFunction16 (RTCOccludedFunc16 occluded16, bool ispc);
      virtual void build(size_t threadIndex, size_t threadCount) {}

    public:
      void* ispcPtr;
      void* ispcIntersect1;
      void* ispcIntersect4;
      void* ispcIntersect8;
      void* ispcIntersect16;
      void* ispcOccluded1;
      void* ispcOccluded4;
      void* ispcOccluded8;
      void* ispcOccluded16;
    };

    struct Instance : public Base
    {
    public:
      Instance (Scene* parent, Accel* object); 
      virtual void setTransform(AffineSpace3f& local2world);
      virtual void build(size_t threadIndex, size_t threadCount) {}
    
    public:
      AffineSpace3f local2world;
      AffineSpace3f world2local;
      Accel* object;
    };
  }
}

#endif
