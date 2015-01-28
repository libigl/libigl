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

#ifndef __EMBREE_GEOMETRY_H__
#define __EMBREE_GEOMETRY_H__

#include "embree2/rtcore.h"
#include "common/default.h"

namespace embree
{
  class Scene;

  /*! type of geometry */
  enum GeometryTy { TRIANGLE_MESH, USER_GEOMETRY, QUADRATIC_BEZIER_CURVES, INSTANCES };
  
#if defined(__SSE__)
  typedef void (*ISPCFilterFunc4)(void* ptr, RTCRay4& ray, __m128 valid);
#endif

#if defined(__AVX__)
  typedef void (*ISPCFilterFunc8)(void* ptr, RTCRay8& ray, __m256 valid);
#endif

#if defined(__MIC__)
  typedef void (*ISPCFilterFunc16)(void* ptr, RTCRay16& ray, __mmask16 valid);
#endif

  /*! Base class all geometries are derived from */
  class Geometry
  {
  public:

    /*! state of a mesh */
    enum State { ENABLING, ENABLED, MODIFIED, DISABLING, DISABLED, ERASING };

  public:
    
    /*! Geometry constructor */
    Geometry (Scene* scene, GeometryTy type, size_t numPrimitives, RTCGeometryFlags flags);

    /*! Virtual destructor */
    virtual ~Geometry() {}

  public:
    __forceinline bool isEnabled() const { 
      return numPrimitives && (state >= ENABLING) && (state <= MODIFIED); 
    }

    __forceinline bool isModified() const { 
      return numPrimitives && ((state == MODIFIED) || (state == ENABLING)); 
    }

    /* test if this is a static geometry */
    __forceinline bool isStatic() const { return flags == RTC_GEOMETRY_STATIC; }

    /* test if this is a deformable geometry */
    __forceinline bool isDeformable() const { return flags == RTC_GEOMETRY_DEFORMABLE; }

    /* test if this is a dynamic geometry */
    __forceinline bool isDynamic() const { return flags == RTC_GEOMETRY_DYNAMIC; }

    /*! for all geometries */
  public:

    /*! Enable geometry. */
    virtual void enable ();

    /*! Update geometry. */
    virtual void update ();
    
    /*! Disable geometry. */
    virtual void disable ();

    /*! Deletes the geometry again. */
    virtual void erase ();

    /*! Free buffers that are unused */
    virtual void immutable () {}

    /*! Verify the geometry */
    virtual bool verify () { return true; }

    /*! called if geometry is switching from disabled to enabled state */
    virtual void enabling() = 0;

    /*! called if geometry is switching from enabled to disabled state */
    virtual void disabling() = 0;

    /*! for triangle mesh and bezier curves only */
  public:

    /*! Sets ray mask. */
    virtual void setMask (unsigned mask) { 
      recordError(RTC_INVALID_OPERATION); 
    }

    /*! Maps specified buffer. */
    virtual void* map(RTCBufferType type) { 
      recordError(RTC_INVALID_OPERATION); 
      return NULL; 
    }

    /*! Unmap specified buffer. */
    virtual void unmap(RTCBufferType type) { 
      recordError(RTC_INVALID_OPERATION); 
    }

    /*! Sets specified buffer. */
    virtual void setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride) { 
      recordError(RTC_INVALID_OPERATION); 
    }

    /*! Set intersection filter function for single rays. */
    virtual void setIntersectionFilterFunction (RTCFilterFunc filter, bool ispc = false);
    
    /*! Set intersection filter function for ray packets of size 4. */
    virtual void setIntersectionFilterFunction4 (RTCFilterFunc4 filter4, bool ispc = false);
    
    /*! Set intersection filter function for ray packets of size 8. */
    virtual void setIntersectionFilterFunction8 (RTCFilterFunc8 filter8, bool ispc = false);
    
    /*! Set intersection filter function for ray packets of size 16. */
    virtual void setIntersectionFilterFunction16 (RTCFilterFunc16 filter16, bool ispc = false);

    /*! Set occlusion filter function for single rays. */
    virtual void setOcclusionFilterFunction (RTCFilterFunc filter, bool ispc = false);
    
    /*! Set occlusion filter function for ray packets of size 4. */
    virtual void setOcclusionFilterFunction4 (RTCFilterFunc4 filter4, bool ispc = false);
    
    /*! Set occlusion filter function for ray packets of size 8. */
    virtual void setOcclusionFilterFunction8 (RTCFilterFunc8 filter8, bool ispc = false);
    
    /*! Set occlusion filter function for ray packets of size 16. */
    virtual void setOcclusionFilterFunction16 (RTCFilterFunc16 filter16, bool ispc = false);

    /*! instances only */
  public:
    
    /*! Sets transformation of the instance */
    virtual void setTransform(AffineSpace3f& transform) {
      recordError(RTC_INVALID_OPERATION); 
    };

    /*! user geometry only */
  public:

    /*! Set user data for intersect and occluded functions. */
    virtual void setUserData (void* ptr, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }

    /*! Set bounds function. */
    virtual void setBoundsFunction (RTCBoundsFunc bounds) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set intersect function for single rays. */
    virtual void setIntersectFunction (RTCIntersectFunc intersect, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set intersect function for ray packets of size 4. */
    virtual void setIntersectFunction4 (RTCIntersectFunc4 intersect4, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set intersect function for ray packets of size 8. */
    virtual void setIntersectFunction8 (RTCIntersectFunc8 intersect8, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set intersect function for ray packets of size 16. */
    virtual void setIntersectFunction16 (RTCIntersectFunc16 intersect16, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }

    /*! Set occlusion function for single rays. */
    virtual void setOccludedFunction (RTCOccludedFunc occluded, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set occlusion function for ray packets of size 4. */
    virtual void setOccludedFunction4 (RTCOccludedFunc4 occluded4, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set occlusion function for ray packets of size 8. */
    virtual void setOccludedFunction8 (RTCOccludedFunc8 occluded8, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }
    
    /*! Set occlusion function for ray packets of size 16. */
    virtual void setOccludedFunction16 (RTCOccludedFunc16 occluded16, bool ispc = false) { 
      recordError(RTC_INVALID_OPERATION); 
    }

  public:
    Scene* parent;   //!< pointer to scene this mesh belongs to
    GeometryTy type;
    ssize_t numPrimitives;    //!< number of primitives of this geometry
    unsigned id;       //!< internal geometry ID
    RTCGeometryFlags flags;    //!< flags of geometry
    State state;       //!< state of the geometry 
    void* userPtr;     //!< user pointer

  public:
    RTCFilterFunc intersectionFilter1;
    RTCFilterFunc occlusionFilter1;

    RTCFilterFunc4 intersectionFilter4;
    RTCFilterFunc4 occlusionFilter4;
    void* ispcIntersectionFilter4; // FIXME: this ISPC mode can be encoded more compactly
    void* ispcOcclusionFilter4;

    RTCFilterFunc8 intersectionFilter8;
    RTCFilterFunc8 occlusionFilter8;
    void* ispcIntersectionFilter8;
    void* ispcOcclusionFilter8;

    RTCFilterFunc16 intersectionFilter16;
    RTCFilterFunc16 occlusionFilter16;
    void* ispcIntersectionFilter16;
    void* ispcOcclusionFilter16;

    __forceinline bool hasIntersectionFilter1() const { return intersectionFilter1 != NULL; }
    __forceinline bool hasIntersectionFilter4() const { return intersectionFilter4 != NULL; }
    __forceinline bool hasIntersectionFilter8() const { return intersectionFilter8 != NULL; }
    __forceinline bool hasIntersectionFilter16() const { return intersectionFilter16 != NULL; }

    __forceinline bool hasOcclusionFilter1() const { return occlusionFilter1 != NULL; }
    __forceinline bool hasOcclusionFilter4() const { return occlusionFilter4 != NULL; }
    __forceinline bool hasOcclusionFilter8() const { return occlusionFilter8 != NULL; }
    __forceinline bool hasOcclusionFilter16() const { return occlusionFilter16 != NULL; }
  };
}

#endif
