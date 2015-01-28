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

#ifndef __EMBREE_ACCELSET_H__
#define __EMBREE_ACCELSET_H__

#include "common/default.h"
#include "common/builder.h"

namespace embree
{
  /*! Base class for set of acceleration structures. */
  class AccelSet : public RefCount 
  {
    ALIGNED_CLASS;
  public:

    typedef RTCIntersectFunc IntersectFunc;
    typedef RTCIntersectFunc4 IntersectFunc4;
    typedef RTCIntersectFunc8 IntersectFunc8;
    typedef RTCIntersectFunc16 IntersectFunc16;
    
    typedef RTCOccludedFunc OccludedFunc;
    typedef RTCOccludedFunc4 OccludedFunc4;
    typedef RTCOccludedFunc8 OccludedFunc8;
    typedef RTCOccludedFunc16 OccludedFunc16;

    struct Intersector1
    {
      Intersector1 (ErrorFunc error = NULL) 
      : intersect((IntersectFunc)error), occluded((OccludedFunc)error), name(NULL) {}

      Intersector1 (IntersectFunc intersect, OccludedFunc occluded, const char* name)
      : intersect(intersect), occluded(occluded), name(name) {}
      
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        IntersectFunc intersect;
        OccludedFunc occluded;  
      };
      
      struct Intersector4 
      {
        Intersector4 (ErrorFunc error = NULL) 
        : intersect((IntersectFunc4)error), occluded((OccludedFunc4)error), name(NULL) {}

        Intersector4 (IntersectFunc4 intersect, OccludedFunc4 occluded, const char* name)
        : intersect(intersect), occluded(occluded), name(name) {}
        
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        IntersectFunc4 intersect;
        OccludedFunc4 occluded;
      };
      
      struct Intersector8 
      {
        Intersector8 (ErrorFunc error = NULL) 
        : intersect((IntersectFunc8)error), occluded((OccludedFunc8)error), name(NULL) {}

        Intersector8 (IntersectFunc8 intersect, OccludedFunc8 occluded, const char* name)
        : intersect(intersect), occluded(occluded), name(name) {}
        
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        IntersectFunc8 intersect;
        OccludedFunc8 occluded;
      };
      
      struct Intersector16 
      {
        Intersector16 (ErrorFunc error = NULL) 
        : intersect((IntersectFunc16)error), occluded((OccludedFunc16)error), name(NULL) {}

        Intersector16 (IntersectFunc16 intersect, OccludedFunc16 occluded, const char* name)
        : intersect(intersect), occluded(occluded), name(name) {}
        
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        IntersectFunc16 intersect;
        OccludedFunc16 occluded;
      };
      
    public:
      
      /*! Construction */
      AccelSet (size_t numItems) : numItems(numItems) {
        intersectors.ptr = NULL; 
        intersectors.boundsPtr = NULL;
      }
      
      /*! Virtual destructor */
      virtual ~AccelSet() {}
      
      /*! makes the acceleration structure immutable */
      virtual void immutable () {};
      
      /*! build accel */
      virtual void build (size_t threadIndex, size_t threadCount) = 0;

      /*! return number of items in set */
      __forceinline size_t size() const {
        return numItems;
      }

      /*! Calculates the bounds of an item */
      __forceinline BBox3f bounds (size_t item) 
      {
        BBox3f box; 
        boundsFunc(intersectors.boundsPtr,item,(RTCBounds&)box);
        return box;
      }
      
      /*! Intersects a single ray with the scene. */
      __forceinline void intersect (RTCRay& ray, size_t item) {
        assert(intersectors.intersector1.intersect);
        intersectors.intersector1.intersect(intersectors.ptr,ray,item);
      }
      
      /*! Intersects a packet of 4 rays with the scene. */
      __forceinline void intersect4 (const void* valid, RTCRay4& ray, size_t item) {
        assert(intersectors.intersector4.intersect);
        intersectors.intersector4.intersect(valid,intersectors.ptr,ray,item);
      }
      
      /*! Intersects a packet of 8 rays with the scene. */
      __forceinline void intersect8 (const void* valid, RTCRay8& ray, size_t item) {
        assert(intersectors.intersector8.intersect);
        intersectors.intersector8.intersect(valid,intersectors.ptr,ray,item);
      }

      /*! Intersects a packet of 16 rays with the scene. */
      __forceinline void intersect16 (const void* valid, RTCRay16& ray, size_t item) {
        assert(intersectors.intersector16.intersect);
        intersectors.intersector16.intersect(valid,intersectors.ptr,ray,item);
      }
      
      /*! Tests if single ray is occluded by the scene. */
      __forceinline void occluded (RTCRay& ray, size_t item) {
        assert(intersectors.intersector1.occluded);
        intersectors.intersector1.occluded(intersectors.ptr,ray,item);
      }
      
      /*! Tests if a packet of 4 rays is occluded by the scene. */
      __forceinline void occluded4 (const void* valid, RTCRay4& ray, size_t item) {
        assert(intersectors.intersector4.occluded);
        intersectors.intersector4.occluded(valid,intersectors.ptr,ray,item);
      }
      
      /*! Tests if a packet of 8 rays is occluded by the scene. */
      __forceinline void occluded8 (const void* valid, RTCRay8& ray, size_t item) {
        assert(intersectors.intersector8.occluded);
        intersectors.intersector8.occluded(valid,intersectors.ptr,ray,item);
      }
      
      /*! Tests if a packet of 16 rays is occluded by the scene. */
      __forceinline void occluded16 (const void* valid, RTCRay16& ray, size_t item) {
        assert(intersectors.intersector16.occluded);
        intersectors.intersector16.occluded(valid,intersectors.ptr,ray,item);
      }
      
    public:
      size_t numItems;
      RTCBoundsFunc boundsFunc;

      struct Intersectors 
      {
        Intersectors() 
          : ptr(NULL), boundsPtr(NULL) {}
      public:
        void* ptr;
        void* boundsPtr;
        Intersector1 intersector1;
        Intersector4 intersector4;
        Intersector8 intersector8;
        Intersector16 intersector16;
      } intersectors;
  };


  struct AccelSetItem {
    AccelSet* accel;
    size_t item;
  };


#define DEFINE_SET_INTERSECTOR1(symbol,intersector)                     \
  AccelSet::Intersector1 symbol((AccelSet::IntersectFunc)intersector::intersect, \
                                (AccelSet::OccludedFunc )intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol));

#define DEFINE_SET_INTERSECTOR4(symbol,intersector)                         \
  AccelSet::Intersector4 symbol((AccelSet::IntersectFunc4)intersector::intersect, \
                                (AccelSet::OccludedFunc4)intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol));

#define DEFINE_SET_INTERSECTOR8(symbol,intersector)                         \
  AccelSet::Intersector8 symbol((AccelSet::IntersectFunc8)intersector::intersect, \
                                (AccelSet::OccludedFunc8)intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol));

#define DEFINE_SET_INTERSECTOR16(symbol,intersector)                         \
  AccelSet::Intersector16 symbol((AccelSet::IntersectFunc16)intersector::intersect, \
                                 (AccelSet::OccludedFunc16)intersector::occluded, \
                                 TOSTRING(isa) "::" TOSTRING(symbol));  
}

#endif
