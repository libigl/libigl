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

#include "sys/platform.h"
#include "sys/taskscheduler.h"
#include "embree/include/embree.h"

#define CATCH_BEGIN try {
#define CATCH_END } catch (std::runtime_error& e) {     \
    std::cerr << "Error: " << e.what() << std::endl;            \
    exit(1);                                                    \
  } catch (...) {                                               \
    std::cerr << "Unknown exception caught." << std::endl;      \
    exit(1);                                                    \
  }

namespace embree
{
  extern "C" void rtcISPCInit() {
    CATCH_BEGIN;
    rtcInit();
    CATCH_END;
  }

  extern "C" void rtcISPCStartThreads(int numThreads) {
    CATCH_BEGIN;
    rtcStartThreads(numThreads);
    CATCH_END;
  }

 extern "C" void rtcISPCStopThreads() {
    CATCH_BEGIN;
    rtcStopThreads();
    CATCH_END;
  }

  extern "C" void rtcISPCExit() {
    CATCH_BEGIN;
    rtcExit();
    CATCH_END;
  }

  extern "C" void rtcISPCDebug() {
    CATCH_BEGIN;
    rtcDebug();
    CATCH_END;
  }

  extern "C" void rtcISPCFreeMemory() {
    CATCH_BEGIN;
    rtcFreeMemory();
    CATCH_END;
  }

  extern "C" void rtcISPCSetVerbose(int verbose) {
    CATCH_BEGIN;
    rtcSetVerbose(verbose);
    CATCH_END;
  }

  extern "C" void* rtcISPCNewVirtualGeometry (const int numObjects) {
    CATCH_BEGIN;
    return rtcNewVirtualGeometry(numObjects,"default");
    CATCH_END;
  }

  extern "C" void rtcISPCSetVirtualGeometryUserData (RTCGeometry* geom, const int i, const int id0, const int id1, const int mask) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryUserData((RTCGeometry*)geom,i,id0,id1,mask);
    CATCH_END;
  }

  extern "C" void rtcISPCSetVirtualGeometryBounds (RTCGeometry* geom, const int i, const void* lower, const void* upper, const void* local2world) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryBounds((RTCGeometry*)geom,i,(const float*)lower,(const float*)upper,(const RTCTransformation*)local2world);
    CATCH_END;
  }

  extern "C" void rtcISPCSetVirtualGeometryIntersector1 (RTCGeometry* geom, const int i, void* intersector1) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryIntersector1((RTCGeometry*)geom,i,(Intersector1*)intersector1);
    CATCH_END;
  }

#if defined (__SSE__)
  extern "C" void rtcISPCSetVirtualGeometryIntersector4 (RTCGeometry* geom, const int i, void* intersector4) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryIntersector4((RTCGeometry*)geom,i,(Intersector4*)intersector4);
    CATCH_END;
  }
#endif

#if defined (__AVX__)
  extern "C" void rtcISPCSetVirtualGeometryIntersector8 (RTCGeometry* geom, const int i, void* intersector8) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryIntersector8((RTCGeometry*)geom,i,(Intersector8*)intersector8);
    CATCH_END;
  }
#endif

#if defined (__MIC__)
  extern "C" void rtcISPCSetVirtualGeometryIntersector16 (RTCGeometry* geom, const int i, void* intersector16) {
    CATCH_BEGIN;
    rtcSetVirtualGeometryIntersector16((RTCGeometry*)geom,i,(Intersector16*)intersector16);
    CATCH_END;
  }
#endif

  extern "C" void* rtcISPCNewTriangleMesh (const int numTriangles, const int numPositions) {
    CATCH_BEGIN;
    return rtcNewTriangleMesh(numTriangles, numPositions, "default");
    CATCH_END;
  }

  extern "C" void* rtcISPCMapPositionBuffer(void* mesh) {
    CATCH_BEGIN;
    return rtcMapPositionBuffer((RTCGeometry*)mesh);
    CATCH_END;
  }

  extern "C" void rtcISPCUnmapPositionBuffer(void* mesh) {
    CATCH_BEGIN;
    rtcUnmapPositionBuffer((RTCGeometry*)mesh);
    CATCH_END;
  }

  extern "C" void* rtcISPCMapTriangleBuffer(void* mesh) {
    CATCH_BEGIN;
    return rtcMapTriangleBuffer((RTCGeometry*)mesh);
    CATCH_END;
  }

  extern "C" void rtcISPCUnmapTriangleBuffer(void* mesh) {
    CATCH_BEGIN;
    rtcUnmapTriangleBuffer((RTCGeometry*)mesh);
    CATCH_END;
  }

  extern "C" void rtcISPCSetApproxBounds (void* geom, const float* lower, const float* upper) {
    CATCH_BEGIN;
    rtcSetApproxBounds((RTCGeometry*)geom,lower,upper);
    CATCH_END;
  }

  extern "C" void rtcISPCBuildAccel (int threadIndex, void* geom) {
    CATCH_BEGIN;
    TaskScheduler::Event* event = TaskScheduler::getISPCEvent(threadIndex);
    rtcBuildAccelAsync((RTCEvent*)event,(RTCGeometry*)geom,"default");
    CATCH_END;
  }

  extern "C" void rtcISPCCleanupGeometry (void* mesh) {
    CATCH_BEGIN;
    rtcCleanupGeometry((RTCGeometry*)mesh);
    CATCH_END;
  }

  extern "C" void rtcISPCDeleteGeometry (void* mesh) {
    CATCH_BEGIN;
    rtcDeleteGeometry((RTCGeometry*)mesh);
    CATCH_END;
  }
  
  extern "C" void rtcISPCGetBounds (void* geom, void* lower_o, void* upper_o) {
    CATCH_BEGIN;
    rtcGetBounds((RTCGeometry*)geom,(float*)lower_o,(float*)upper_o);
    CATCH_END;
  }

  extern "C" void* rtcISPCQueryIntersector1 (const void* geom, const void* type) {
    CATCH_BEGIN;
    if (type == NULL) return rtcQueryIntersector1((const RTCGeometry*)geom,"default");
    else return rtcQueryIntersector1((const RTCGeometry*)geom,(const char*)type);
    CATCH_END;
  }

  extern "C" void rtcISPCDeleteIntersector1 (void* intersector) {
    CATCH_BEGIN;
    rtcDeleteIntersector1((Intersector1*)intersector);
    CATCH_END;
  }

#if defined (__SSE__)

  extern "C" void* rtcISPCQueryIntersector4 (const void* geom, const void* type) 
  {
    CATCH_BEGIN;
    if (type == NULL) return rtcQueryIntersector4((const RTCGeometry*)geom,"default");
    else              return rtcQueryIntersector4((const RTCGeometry*)geom,(const char*)type);
    CATCH_END;
  }

  extern "C" void rtcISPCDeleteIntersector4 (void* intersector) 
  {
    CATCH_BEGIN;
    rtcDeleteIntersector4((Intersector4*)intersector);
    CATCH_END;
  }

#endif

#if defined (__AVX__)

  extern "C" void* rtcISPCQueryIntersector8 (const void* geom, const void* type) 
  {
    CATCH_BEGIN;
    if (type == NULL) return rtcQueryIntersector8((const RTCGeometry*)geom,"default");
    else              return rtcQueryIntersector8((const RTCGeometry*)geom,(const char*)type);
    CATCH_END;
  }

  extern "C" void rtcISPCDeleteIntersector8 (void* intersector) 
  {
    CATCH_BEGIN;
    rtcDeleteIntersector8((Intersector8*)intersector);
    CATCH_END;
  }

#endif

#if defined (__MIC__)

  extern "C" void* rtcISPCQueryIntersector16 (const void* geom, const void* type) 
  {
    CATCH_BEGIN;
    if (type == NULL) return rtcQueryIntersector16((const RTCGeometry*)geom,"default");
    else              return rtcQueryIntersector16((const RTCGeometry*)geom,(const char*)type);
    CATCH_END;
  }

  extern "C" void rtcISPCDeleteIntersector16 (void* intersector) 
  {
    CATCH_BEGIN;
    rtcDeleteIntersector16((Intersector16*)intersector);
    CATCH_END;
  }

#endif
}
