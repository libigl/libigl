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

#include "include/embree.h"
#include "common/alloc.h"
#include "common/registry_intersector.h"
#include "common/registry_builder.h"
#include "geometry/triangle_mesh.h"
#include "geometry/virtual_scene.h"
#include "geometry/triangles.h"

namespace embree
{
  /*! verbose mode */
  int g_verbose = 0;

  /*! register functions for acceleration structures */
  void BVH2Register();
  void BVH4Register();
  void BVH4MBRegister();
  void BVH4iRegister();
  void BVH8Register();
  void BVH4AOSRegister();

  /*! primitive types */
  const VirtualObject ::Type VirtualObject ::type;
  const Triangle1i::Type Triangle1i::type;
  const Triangle1v::Type Triangle1v::type;
  const Triangle1 ::Type Triangle1 ::type;

#if defined (__SSE__)
  const Triangle4i::Type Triangle4i::type;
  const Triangle4v::Type Triangle4v::type;
  const Triangle4 ::Type Triangle4 ::type;
#endif

#if defined (__AVX__)
  const Triangle8 ::Type Triangle8 ::type;
#endif

  /*! intersector interface names */
  const char* const Intersector1::name = "Intersector1";

#if defined(__SSE__)
  const char* const Intersector4::name = "Intersector4";
#endif

#if defined(__AVX__)
  const char* const Intersector8::name = "Intersector8";
#endif

#if defined(__MIC__)
  const char* const Intersector16::name = "Intersector16";
#endif

  static bool initialized = false;
  
  void rtcInit() 
  {
    if (initialized) 
      throw std::runtime_error("embree is already initialized");

    static bool registered = false;
    if (!registered) 
    {
      TriangleMesh::accels.add("triangle1i",Triangle1i::type);
      TriangleMesh::accels.add("triangle1v",Triangle1v::type);
      TriangleMesh::accels.add("triangle1" ,Triangle1 ::type);

#if defined (__SSE__)
      TriangleMesh::accels.add("triangle4i",Triangle4i::type);
      TriangleMesh::accels.add("triangle4v",Triangle4v::type);
      TriangleMesh::accels.add("triangle4" ,Triangle4 ::type);
#endif

#if defined (__AVX__)
      TriangleMesh::accels.add("triangle8",Triangle8::type);
#endif

      VirtualScene::accels.add("virtual",VirtualObject::type);

      BVH2Register();
      
#if defined (__SSE__)
      BVH4Register();
      BVH4MBRegister();
      BVH4iRegister();
#endif
      
#if defined (__AVX__)
      BVH8Register();
#endif
      
#if defined (__MIC__)
      BVH4AOSRegister();
#endif
    }

    initialized = true;
  }

  void rtcExit() 
  {
    if (!initialized)
      throw std::runtime_error("embree is not initialized");

    TriangleMesh::clearRegistry();
    VirtualScene::clearRegistry();
    TaskScheduler::destroy();
    initialized = false;
  }

  void rtcDebug() 
  {
#if defined(__USE_STAT_COUNTERS__)
    Stat::print(std::cout);
    Stat::clear();
#endif
  }

  void rtcStartThreads(size_t numThreads) {
    TaskScheduler::create(numThreads);
  }

  void rtcStopThreads() {
    TaskScheduler::destroy();
  }

  void rtcFreeMemory() {
    Alloc::global.clear();
  }

  void rtcSetVerbose(int verbose) {
    g_verbose = verbose;
  }

  RTCGeometry* rtcNewTriangleMesh (const size_t numTriangles, const size_t numPositions, const char* accelTy) {
    return new TriangleMesh(numTriangles,numPositions,empty,accelTy);
  }

  RTCVertex* rtcMapPositionBuffer(RTCGeometry* mesh) {
    return ((TriangleMesh*)mesh)->vertices;
  }

  void rtcUnmapPositionBuffer(RTCGeometry* mesh) {
  }

  RTCTriangle* rtcMapTriangleBuffer(RTCGeometry* mesh) {
    return ((TriangleMesh*)mesh)->triangles;
  }

  void rtcUnmapTriangleBuffer(RTCGeometry* mesh) {
  }

  void rtcCleanupGeometry (RTCGeometry* mesh) {
    mesh->freeze();
  }

  RTCGeometry* rtcNewVirtualGeometry (const size_t numObjects, const char* accelTy) {
    return new VirtualScene(numObjects,accelTy);
  }

  void rtcSetVirtualGeometryUserData (RTCGeometry* geom, const size_t i, const int id0, const int id1, const int mask)
  {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    obj.id0 = id0;
    obj.id1 = id1;
    obj.mask = mask;
  }

  void rtcSetVirtualGeometryBounds (RTCGeometry* geom, const size_t i, const float* lower, const float* upper, const RTCTransformation* local2world)
  {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    if (lower && upper) {
      obj.localBounds = BBox3f(Vector3f(lower[0],lower[1],lower[2]),
                               Vector3f(upper[0],upper[1],upper[2]));
    }
    if (local2world) {
      obj.local2world = AffineSpace3f(Vector3f(local2world->vxx,local2world->vxy,local2world->vxz),
                                      Vector3f(local2world->vyx,local2world->vyy,local2world->vyz),
                                      Vector3f(local2world->vzx,local2world->vzy,local2world->vzz),
                                      Vector3f(local2world->px ,local2world->py ,local2world->pz ));
      obj.hasTransform = obj.local2world != AffineSpace3f(one);
    }
    obj.calculateWorldData();
  }

  void rtcSetVirtualGeometryIntersector1 (RTCGeometry* geom, const size_t i, const Intersector1* intersector1) {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    obj.intersector1 = intersector1;
  }

#if defined(__SSE__)
  void rtcSetVirtualGeometryIntersector4 (RTCGeometry* geom, const size_t i, const Intersector4* intersector4) {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    obj.intersector4 = intersector4;
  }
#endif

#if defined(__AVX__)
  void rtcSetVirtualGeometryIntersector8 (RTCGeometry* geom, const size_t i, const Intersector8* intersector8) {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    obj.intersector8 = intersector8;
  }
#endif

#if defined(__MIC__)
  void rtcSetVirtualGeometryIntersector16 (RTCGeometry* geom, const size_t i, const Intersector16* intersector16) {
    VirtualScene::Object& obj = ((VirtualScene*)geom)->get(i);
    obj.intersector16 = intersector16;
  }
#endif

  void rtcSetApproxBounds (RTCGeometry* geom, const float* lower, const float* upper) 
  {
    geom->approx.lower = Vector3f(lower[0],lower[1],lower[2]);
    geom->approx.upper = Vector3f(upper[0],upper[1],upper[2]);
  }

  void rtcDeleteGeometry(RTCGeometry* mesh) {
    delete mesh;
  }

  void rtcBuildAccel(RTCGeometry* geom, const char* builderTy)
  {
    double t0 = getSeconds();
    TaskScheduler::EventSync event;
    geom->build(&event,builderTy);
    event.sync();
    double dt = getSeconds()-t0;

    /* output statistics */
    if (g_verbose > 0) {
      std::ostringstream stream;
      size_t numTriangles = geom->size();
      stream.setf(std::ios::fixed, std::ios::floatfield);
      stream.precision(0);
      stream << "triangles = " << numTriangles << std::endl;
      stream << "build time = " << dt*1000.0f << " ms" << std::endl;
      stream.precision(3);
      stream << "build performance = " << double(numTriangles)/dt*1E-6 << " Mtris/s" << std::endl;
      stream << "memory pool = " << double(Alloc::global.size())*1E-6 << " MB" << std::endl;
      std::cout << stream.str();
      geom->accel->print();
    }
  }

  void rtcBuildAccelAsync (RTCEvent* event, RTCGeometry* geom, const char* builderTy) {
    geom->build((TaskScheduler::Event*)event,builderTy);
  }

  void rtcGetBounds (RTCGeometry* geom, float* lower_o, float* upper_o)
  {
    BBox3f bounds = geom->bounds();
    lower_o[0] = bounds.lower.x; lower_o[1] = bounds.lower.y; lower_o[2] = bounds.lower.z;
    upper_o[0] = bounds.upper.x; upper_o[1] = bounds.upper.y; upper_o[2] = bounds.upper.z;
  }

  Intersector1* rtcQueryIntersector1 (const RTCGeometry* geom, const char* type) {
    return geom->intersector1(type);
  }

  void rtcDeleteIntersector1 (Intersector1* intersector) {
    delete intersector;
  }
 
#if defined(__SSE__)

  Intersector4* rtcQueryIntersector4 (const RTCGeometry* geom, const char* type) {
    return geom->intersector4(type);
  }

  void rtcDeleteIntersector4 (Intersector4* intersector) {
    delete intersector;
  }

#endif

#if defined(__AVX__)

  Intersector8* rtcQueryIntersector8 (const RTCGeometry* geom, const char* type) {
    return geom->intersector8(type);
  }

  void rtcDeleteIntersector8 (Intersector8* intersector) {
    delete intersector;
  }

#endif

#if defined(__MIC__)

  Intersector16* rtcQueryIntersector16 (const RTCGeometry* geom, const char* type) {
    return geom->intersector16(type);
  }

  void rtcDeleteIntersector16 (Intersector16* intersector) {
    delete intersector;
  }

#endif

  RTCEvent* rtcNewEvent() {
    return (RTCEvent*) new TaskScheduler::EventSync;
  }

  void rtcWait(RTCEvent* event) {
    ((TaskScheduler::EventSync*) event)->sync();
  }

  void rtcDeleteEvent(RTCEvent* event) {
    delete ((TaskScheduler::EventSync*) event);
  }
}

