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

#ifndef __EMBREE_RTCORE_H__
#define __EMBREE_RTCORE_H__

#if defined(_MSC_VER) && !defined(__SSE__)
#define __SSE__
#endif

#if !defined(_WIN32) && !defined(__forceinline)
#define __forceinline inline __attribute__((always_inline))
#endif

#include "intersector1.h"

#if defined(__SSE__)
#include "intersector4.h"
#endif

#if defined(__AVX__)
#include "intersector8.h"
#endif

#if defined(__MIC__)
#include "intersector16.h"
#endif

namespace embree
{
  /*! Opaque event type */
  struct RTCEvent;

  /*! Opaque geometry type */
  struct RTCGeometry;

  /*! Triangle indices */
  struct RTCTriangle
  {
    /*! Constructs a builder triangle. */
    __forceinline RTCTriangle(int v0,      //!< reference to 1st vertex of the triangle
                              int v1,      //!< reference to 2nd vertex of the triangle
                              int v2,      //!< reference to 3rd vertex of the triangle
                              int id0 = 0, //!< 1st optional user ID
                              int id1 = 0) //!< 2nd optional user ID
      : v0(v0), v1(v1), v2(v2), id0(id0), id1(id1) { }
    
  public:
    int v0, v1, v2;       //!< references to 1st, 2nd, and 3rd triangle vertex
    int id0, id1;         //!< two user IDs
  };
  
  /*! Vertex type */
  struct RTCVertex {
  public:
    RTCVertex (float x, float y, float z) : x(x), y(y), z(z), align(0) {}
  public:
    float x,y,z; 
  private:
    int align;
  };

  /*! Transformation type */
  struct RTCTransformation 
  {
    float vxx,vxy,vxz; 
    float vyx,vyy,vyz; 
    float vzx,vzy,vzz; 
    float px ,py ,pz ; 
  };
  
  /*! Initialized the ray tracing core. */
  void rtcInit();

  /*! Starts the Embree threads. Embree threads might spin depending on implementation. */
  void rtcStartThreads(size_t numThreads = 0);

  /*! Stops the Embree threads. Threads will stop spinning. */
  void rtcStopThreads();

  /*! Cleanup the ray tracing core. */
  void rtcExit();

  /*! Print debug info. */
  void rtcDebug();

  /*! Frees all unused internal memory. */
  void rtcFreeMemory();

  /*! enable verbose output */
  void rtcSetVerbose(int verbose);

  /*! New list of objects. */
  RTCGeometry* rtcNewVirtualGeometry (const size_t numObjects,  //!< maximal number of objects
                                      const char* accelTy       //!< type of acceleration structure to use
                                      );

  /*! Set virtual object user data. */
  void rtcSetVirtualGeometryUserData (RTCGeometry* geom, const size_t i, const int id0, const int id1, const int mask = -1);

  /*! Set virtual object bounding box and transformation. */
  void rtcSetVirtualGeometryBounds (RTCGeometry* geom, const size_t i, const float* lower, const float* upper, const RTCTransformation* local2world = NULL);

  /*! Set virtual object intersector1. */
  void rtcSetVirtualGeometryIntersector1 (RTCGeometry* geom, const size_t i, const RTCIntersector1* intersector1);

  /*! Set virtual object intersector4. */
#if defined(__SSE__)
  void rtcSetVirtualGeometryIntersector4 (RTCGeometry* geom, const size_t i, const RTCIntersector4* intersector4);
#endif

  /*! Set virtual object intersector8. */
#if defined(__AVX__)
  void rtcSetVirtualGeometryIntersector8 (RTCGeometry* geom, const size_t i, const RTCIntersector8* intersector8);
#endif

  /*! Set virtual object intersector16. */
#if defined(__MIC__)
  void rtcSetVirtualGeometryIntersector16 (RTCGeometry* geom, const size_t i, const RTCIntersector16* intersector16);
#endif

  /*! Creates new triangle mesh. */
  RTCGeometry* rtcNewTriangleMesh (const size_t numTriangles, //!< number of triangles
                                   const size_t numPositions, //!< maximal number of vertices
                                   const char* accelTy        //!< type of acceleration structure to use
                                   );

  /*! Map position buffer */
  RTCVertex* rtcMapPositionBuffer(RTCGeometry* mesh); 

  /*! Unmap position buffer */
  void rtcUnmapPositionBuffer(RTCGeometry* mesh); 

  /*! Map triangle buffer. */
  RTCTriangle* rtcMapTriangleBuffer(RTCGeometry* mesh);

  /*! Unmap triangle buffer. */
  void rtcUnmapTriangleBuffer(RTCGeometry* mesh); 

  /*! Set optional approximate bounds. */
  void rtcSetApproxBounds (RTCGeometry* geom, const float* lower, const float* upper);

  /*! Builds acceleration structure of specified type over some geometry. */
  void rtcBuildAccel (RTCGeometry* geom,     //!< geometry to build acceleration structure for
                      const char* builderTy  //!< builder to use
                      ); 

  /*! Asynchronous build of acceleration structure of specified type over some geometry. */
  void rtcBuildAccelAsync (RTCEvent* event,        //!< event to trigger if build is finished
                           RTCGeometry* geom,      //!< geometry to build acceleration structure for
                           const char* builderTy   //!< builder to use
                           ); 

  /*! Free unused data of the geometry. A later modification of the geometry is no longer possible. */
  void rtcCleanupGeometry (RTCGeometry* mesh);

  /*! Destroys geometry. */
  void rtcDeleteGeometry (RTCGeometry* mesh);

  /*! Returns bounding box. */
  void rtcGetBounds (RTCGeometry* geom, float* lower_o, float* upper_o);

  /*! Query RTCIntersector1 interface for single rays */
  RTCIntersector1* rtcQueryIntersector1 (const RTCGeometry* geom, const char* type);

  /*! Delete intersector1. */
  void rtcDeleteIntersector1 (RTCIntersector1* intersector);

#if defined(__SSE__)

  /*! Query Intersector4 interface for packets of 4 rays. */
  RTCIntersector4* rtcQueryIntersector4 (const RTCGeometry* geom, const char* type);
  
  /*! Delete intersector4. */
  void rtcDeleteIntersector4 (RTCIntersector4* intersector);

#endif

#if defined(__AVX__)

  /*! Query RTCIntersector8 interface for packets of 8 rays. */
  RTCIntersector8* rtcQueryIntersector8 (const RTCGeometry* geom, const char* type);

  /*! Delete intersector8. */
  void rtcDeleteIntersector8 (RTCIntersector8* intersector);

#endif

#if defined(__MIC__)

  /*! Query RTCIntersector16 interface for packets of 16 rays. */
  RTCIntersector16* rtcQueryIntersector16 (const RTCGeometry* geom, const char* type);

  /*! Delete intersector16. */
  void rtcDeleteIntersector16 (RTCIntersector16* intersector);

#endif

  /*! Create a new event */
  RTCEvent* rtcNewEvent();

  /*! Waits for an event */
  void rtcWait(RTCEvent* event);

  /*! Delete event */
  void rtcDeleteEvent(RTCEvent* event);
}

#endif
