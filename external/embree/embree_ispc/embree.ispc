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

#include "include/embree.isph"

extern "C" void rtcISPCInit();
extern "C" void rtcISPCExit();
extern "C" void rtcISPCDebug();
extern "C" void rtcISPCStartThreads(uniform int numThreads);
extern "C" void rtcISPCStopThreads();
extern "C" void rtcISPCFreeMemory();
extern "C" void rtcISPCSetVerbose(uniform int verbose);

extern "C" void* uniform rtcISPCNewVirtualGeometry (const uniform int numObjects);
extern "C" void rtcISPCSetVirtualGeometryUserData (void* uniform geom, const uniform int i, const uniform int id0, const uniform int id1, const uniform int mask);
extern "C" void rtcISPCSetVirtualGeometryBounds (void* uniform geom, const uniform int i, 
                                                 const uniform float* uniform lower, const uniform float* uniform upper, 
                                                 const void* uniform local2world);
extern "C" void rtcISPCSetVirtualGeometryIntersector4 (void* uniform geom, const uniform int i, void* uniform intersector4);
extern "C" void rtcISPCSetVirtualGeometryIntersector8 (void* uniform geom, const uniform int i, void* uniform intersector8);
extern "C" void rtcISPCSetVirtualGeometryIntersector16 (void* uniform geom, const uniform int i, void* uniform intersector16);

extern "C" void* uniform rtcISPCNewTriangleMesh (const uniform int numTriangles, const uniform int numPositions);
extern "C" void* uniform rtcISPCMapPositionBuffer(void* uniform mesh); 
extern "C" void rtcISPCUnmapPositionBuffer(void* uniform mesh); 
extern "C" void* uniform rtcISPCMapTriangleBuffer(void* uniform mesh);
extern "C" void rtcISPCUnmapTriangleBuffer(void* uniform mesh); 

extern "C" void rtcISPCSetApproxBounds (void* uniform geom, const uniform float* uniform lower, const uniform float* uniform upper);
extern "C" void rtcISPCBuildAccel (uniform int threadIndex, void* uniform geom);
extern "C" void rtcISPCCleanupGeometry (void* uniform mesh);
extern "C" void rtcISPCDeleteGeometry (void* uniform mesh);

extern "C" void rtcISPCGetBounds (void* uniform geom, void* uniform lower_o, void* uniform upper_o);
extern "C" void* uniform rtcISPCQueryIntersector4 (const void* uniform accel, const void* uniform type);
extern "C" void* uniform rtcISPCQueryIntersector8 (const void* uniform accel, const void* uniform type);
extern "C" void* uniform rtcISPCQueryIntersector16(const void* uniform accel, const void* uniform type);
extern "C" void rtcISPCDeleteIntersector4 (void* uniform intersector);
extern "C" void rtcISPCDeleteIntersector8 (void* uniform intersector);
extern "C" void rtcISPCDeleteIntersector16 (void* uniform intersector);

void rtcInit() {
  rtcISPCInit();
}

void rtcExit() {
  rtcISPCExit();
}

void rtcDebug() {
  rtcISPCDebug();
}

void rtcStartThreads(uniform int numThreads) {
  rtcISPCStartThreads(numThreads);
}

void rtcStartThreads() {
  rtcISPCStartThreads(0);
}

void rtcStopThreads() {
  rtcISPCStopThreads();
}

void rtcFreeMemory() {
  rtcISPCFreeMemory();
}

void rtcSetVerbose(uniform int verbose) {
  rtcISPCSetVerbose(verbose);
}

uniform RTCGeometry* uniform rtcNewVirtualGeometry (const uniform int numObjects) {
  return (uniform RTCGeometry* uniform) rtcISPCNewVirtualGeometry(numObjects);
}

void rtcSetVirtualGeometryUserData (uniform RTCGeometry* uniform geom, const uniform int i, const uniform int id0, const uniform int id1)
{
  rtcISPCSetVirtualGeometryUserData (geom, i, id0, id1, -1);
}

void rtcSetVirtualGeometryUserData (uniform RTCGeometry* uniform geom, const uniform int i, const uniform int id0, const uniform int id1, const uniform int mask)
{
  rtcISPCSetVirtualGeometryUserData (geom, i, id0, id1, mask);
}

void rtcSetVirtualGeometryBounds (uniform RTCGeometry* uniform geom, const uniform int i, 
                                  const uniform float* uniform lower, const uniform float* uniform upper, 
                                  const uniform RTCTransformation* uniform local2world)
{
  rtcISPCSetVirtualGeometryBounds (geom, i, lower, upper, local2world);
}

void rtcSetVirtualGeometryIntersector (uniform RTCGeometry* uniform geom, const uniform int i, const uniform Intersector* uniform intersector) {
#if defined (__SSE4__)
  rtcISPCSetVirtualGeometryIntersector4(geom, i, intersector);
#elif defined (__AVX__)
  rtcISPCSetVirtualGeometryIntersector8(geom, i, intersector);
#elif defined(__MIC__) 
  rtcISPCSetVirtualGeometryIntersector16(geom, i, intersector);
#else
#  error "Unsupported SIMD configuration."
#endif
}

uniform RTCGeometry* uniform rtcNewTriangleMesh (const uniform int numTriangles, const uniform int numPositions) {
  return (uniform RTCGeometry* uniform) rtcISPCNewTriangleMesh (numTriangles, numPositions);
}

uniform RTCVertex* uniform rtcMapPositionBuffer(uniform RTCGeometry* uniform mesh) {
  return (uniform RTCVertex* uniform) rtcISPCMapPositionBuffer (mesh);
}

void rtcUnmapPositionBuffer(uniform RTCGeometry* uniform mesh) {
  rtcISPCUnmapPositionBuffer (mesh);
}

uniform RTCTriangle* uniform rtcMapTriangleBuffer(uniform RTCGeometry* uniform mesh) {
  return (uniform RTCTriangle* uniform) rtcISPCMapTriangleBuffer (mesh);
}

void rtcUnmapTriangleBuffer(uniform RTCGeometry* uniform mesh) {
  rtcISPCUnmapTriangleBuffer (mesh);
}

void rtcSetApproxBounds (uniform RTCGeometry* uniform geom, const uniform float* uniform lower, const uniform float* uniform upper) {
  rtcISPCSetApproxBounds (geom, lower, upper);
}

task void rtcBuildAccel (uniform RTCGeometry* uniform geom) {
  rtcISPCBuildAccel (threadIndex,geom);
}

void rtcCleanupGeometry (uniform RTCGeometry* uniform mesh) {
  rtcISPCCleanupGeometry (mesh);
}

void rtcDeleteGeometry (uniform RTCGeometry* uniform mesh) {
  rtcISPCDeleteGeometry (mesh);
}

void rtcGetBounds (uniform RTCGeometry* uniform geom, uniform float* uniform lower_o, uniform float* uniform upper_o) {
  rtcISPCGetBounds (geom, lower_o, upper_o);
}

uniform Intersector* uniform rtcQueryIntersector (uniform RTCGeometry* uniform geom) {
  return rtcQueryIntersector(geom, NULL);
}

uniform Intersector* uniform rtcQueryIntersector (uniform RTCGeometry* uniform geom, void* uniform type)
{
#if defined (__SSE4__)
  return (uniform Intersector* uniform) rtcISPCQueryIntersector4 (geom,type);
#elif defined (__AVX__)
  return (uniform Intersector* uniform) rtcISPCQueryIntersector8 (geom,type);
#elif defined(__MIC__) 
  return (uniform Intersector* uniform) rtcISPCQueryIntersector16 (geom,type);
#else
#  error "Unsupported SIMD configuration."
#endif
}

void rtcDeleteIntersector (uniform Intersector* uniform intersector)
{
#if defined (__SSE4__)
  rtcISPCDeleteIntersector4 (intersector);
#elif defined (__AVX__)
  rtcISPCDeleteIntersector8 (intersector);
#elif defined(__MIC__) 
  rtcISPCDeleteIntersector16 (intersector);
#else
#  error "Unsupported SIMD configuration."
#endif
}

