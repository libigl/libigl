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

#pragma once

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

#if defined(__SSE__)
    typedef void (*ISPCIntersectFunc4)(void* ptr, RTCRay4& ray, size_t item, __m128 valid);
    typedef void (*ISPCOccludedFunc4 )(void* ptr, RTCRay4& ray, size_t item, __m128 valid);
#endif

#if defined(__AVX__)
    typedef void (*ISPCIntersectFunc8)(void* ptr, RTCRay8& ray, size_t item, __m256 valid);
    typedef void (*ISPCOccludedFunc8 )(void* ptr, RTCRay8& ray, size_t item, __m256 valid);
#endif

#if defined(__MIC__)
    typedef void (*ISPCIntersectFunc16)(void* ptr, RTCRay16& ray, size_t item, __mmask16 valid);
    typedef void (*ISPCOccludedFunc16 )(void* ptr, RTCRay16& ray, size_t item, __mmask16 valid);
#endif

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
        : intersect((void*)error), occluded((void*)error), name(NULL), ispc(false) {}

        Intersector4 (void* intersect, void* occluded, const char* name, bool ispc)
        : intersect(intersect), occluded(occluded), name(name), ispc(ispc) {}
	
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        void* intersect;
        void* occluded;
	bool ispc;
      };
      
      struct Intersector8 
      {
        Intersector8 (ErrorFunc error = NULL) 
        : intersect((void*)error), occluded((void*)error), name(NULL), ispc(false) {}

        Intersector8 (void* intersect, void* occluded, const char* name, bool ispc)
        : intersect(intersect), occluded(occluded), name(name), ispc(ispc) {}
        
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        void* intersect;
        void* occluded;
	bool ispc;
      };
      
      struct Intersector16 
      {
        Intersector16 (ErrorFunc error = NULL) 
        : intersect((void*)error), occluded((void*)error), name(NULL), ispc(false) {}

        Intersector16 (void* intersect, void* occluded, const char* name, bool ispc)
        : intersect(intersect), occluded(occluded), name(name), ispc(ispc) {}
        
        operator bool() const { return name; }
        
      public:
        static const char* type;
        const char* name;
        void* intersect;
        void* occluded;
	bool ispc;
      };
      
    public:
      
      /*! Construction */
      AccelSet (size_t numItems) : numItems(numItems) {
        intersectors.ptr = NULL; 
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
      __forceinline BBox3fa bounds (size_t item) const
      {
        BBox3fa box; 
        boundsFunc(intersectors.ptr,item,(RTCBounds&)box);
        return box;
      }
      
      /*! Intersects a single ray with the scene. */
      __forceinline void intersect (RTCRay& ray, size_t item) {
        assert(intersectors.intersector1.intersect);
        intersectors.intersector1.intersect(intersectors.ptr,ray,item);
      }
      
      /*! Intersects a packet of 4 rays with the scene. */
      __forceinline void intersect4 (const void* valid, RTCRay4& ray, size_t item) {
#if defined(__SSE__)
        assert(intersectors.intersector4.intersect);
	if (intersectors.intersector4.ispc) ((ISPCIntersectFunc4)intersectors.intersector4.intersect)(intersectors.ptr,ray,item,*(__m128*)valid);
        else                                ((    IntersectFunc4)intersectors.intersector4.intersect)(valid,intersectors.ptr,ray,item);
#endif
      }
      
      /*! Intersects a packet of 8 rays with the scene. */
      __forceinline void intersect8 (const void* valid, RTCRay8& ray, size_t item) {
#if defined(__AVX__)
        assert(intersectors.intersector8.intersect);
	if (intersectors.intersector8.ispc) ((ISPCIntersectFunc8)intersectors.intersector8.intersect)(intersectors.ptr,ray,item,*(__m256*)valid);
        else                                ((    IntersectFunc8)intersectors.intersector8.intersect)(valid,intersectors.ptr,ray,item);
#endif
      }

      /*! Intersects a packet of 16 rays with the scene. */
      __forceinline void intersect16 (const void* valid, RTCRay16& ray, size_t item) {
#if defined(__MIC__)
        assert(intersectors.intersector16.occluded);
	if (intersectors.intersector16.ispc) {
	  const mic_i maski = *(mic_i*)valid;
	  const __mmask16 mask = maski != mic_i(0);
	  ((ISPCIntersectFunc16)intersectors.intersector16.intersect)(intersectors.ptr,ray,item,mask);
	}
        else
	  ((IntersectFunc16)intersectors.intersector16.intersect)(valid,intersectors.ptr,ray,item);
#endif
      }
      
      /*! Tests if single ray is occluded by the scene. */
      __forceinline void occluded (RTCRay& ray, size_t item) {
        assert(intersectors.intersector1.occluded);
        intersectors.intersector1.occluded(intersectors.ptr,ray,item);
      }
      
      /*! Tests if a packet of 4 rays is occluded by the scene. */
#if defined(__SSE__)
      __forceinline void occluded4 (const void* valid, RTCRay4& ray, size_t item) {
	assert(intersectors.intersector4.occluded);
	if (intersectors.intersector4.ispc) ((ISPCOccludedFunc4)intersectors.intersector4.occluded)(intersectors.ptr,ray,item,*(__m128*)valid);
        else                                ((    OccludedFunc4)intersectors.intersector4.occluded)(valid,intersectors.ptr,ray,item);
      }
#endif
      
      /*! Tests if a packet of 8 rays is occluded by the scene. */
#if defined(__AVX__)
      __forceinline void occluded8 (const void* valid, RTCRay8& ray, size_t item) {
	assert(intersectors.intersector8.occluded);
	if (intersectors.intersector8.ispc) ((ISPCOccludedFunc8)intersectors.intersector8.occluded)(intersectors.ptr,ray,item,*(__m256*)valid);
        else                                ((    OccludedFunc8)intersectors.intersector8.occluded)(valid,intersectors.ptr,ray,item);
      }
#endif
      
      /*! Tests if a packet of 16 rays is occluded by the scene. */
#if defined(__MIC__)
      __forceinline void occluded16 (const void* valid, RTCRay16& ray, size_t item) {
        assert(intersectors.intersector16.occluded);
	if (intersectors.intersector16.ispc) {
	  const mic_i maski = *(mic_i*)valid;
	  const __mmask16 mask = maski != mic_i(0);
	  ((ISPCOccludedFunc16)intersectors.intersector16.occluded)(intersectors.ptr,ray,item,mask);
	}
	else
	  ((OccludedFunc16)intersectors.intersector16.occluded)(valid,intersectors.ptr,ray,item);
      }
#endif
      
    public:
      size_t numItems;
      RTCBoundsFunc boundsFunc;

      struct Intersectors 
      {
        Intersectors() : ptr(NULL) {}
      public:
        void* ptr;
        Intersector1 intersector1;
        Intersector4 intersector4;
        Intersector8 intersector8;
        Intersector16 intersector16;
      } intersectors;
  };

#define DEFINE_SET_INTERSECTOR1(symbol,intersector)                     \
  AccelSet::Intersector1 symbol((AccelSet::IntersectFunc)intersector::intersect, \
                                (AccelSet::OccludedFunc )intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol));

#define DEFINE_SET_INTERSECTOR4(symbol,intersector)                         \
  AccelSet::Intersector4 symbol((void*)intersector::intersect, \
                                (void*)intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol),	\
				false);

#define DEFINE_SET_INTERSECTOR8(symbol,intersector)                         \
  AccelSet::Intersector8 symbol((void*)intersector::intersect, \
                                (void*)intersector::occluded, \
                                TOSTRING(isa) "::" TOSTRING(symbol),	\
				false);

#define DEFINE_SET_INTERSECTOR16(symbol,intersector)                         \
  AccelSet::Intersector16 symbol((void*)intersector::intersect, \
                                 (void*)intersector::occluded, \
                                 TOSTRING(isa) "::" TOSTRING(symbol),\
				 false);  
}
