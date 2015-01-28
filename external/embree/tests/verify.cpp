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
#include "sys/ref.h"
#include "sys/thread.h"
#include "sys/sysinfo.h"
#include "sys/sync/barrier.h"
#include "sys/sync/mutex.h"
#include "sys/sync/condition.h"
#include "math/vec3.h"
#include "math/bbox.h"
#include "embree2/rtcore.h"
#include "embree2/rtcore_ray.h"
#include "../kernels/common/default.h"
#include <vector>

namespace embree
{
#if !defined(__MIC__)
  RTCAlgorithmFlags aflags = (RTCAlgorithmFlags) (RTC_INTERSECT1 | RTC_INTERSECT4 | RTC_INTERSECT8);
#else
  RTCAlgorithmFlags aflags = (RTCAlgorithmFlags) (RTC_INTERSECT1 | RTC_INTERSECT16);
#endif
  /* configuration */
  static std::string g_rtcore = "";
  static size_t testN = 100000;

  /* vertex and triangle layout */
  struct Vertex   { float x,y,z,a; };
  struct Triangle { int v0, v1, v2; };

  std::vector<thread_t> g_threads;

#define AssertNoError() \
  if (rtcGetError() != RTC_NO_ERROR) return false;
#define AssertAnyError() \
  if (rtcGetError() == RTC_NO_ERROR) return false;
#define AssertError(code) \
  if (rtcGetError() != code) return false;
#define POSITIVE(name,test) {                                               \
    printf("%30s ...",name);                                            \
    bool ok = test;                                                     \
    printf(" %s\n",ok ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m");          \
	fflush(stdout);\
  }
#if defined(__EXIT_ON_ERROR__)
#define NEGATIVE(name,test)
#else
#define NEGATIVE(name,test) {                                       \
    printf("%30s ... ",name);                                           \
    bool notok = test;                                                  \
    printf(" %s\n",notok ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m"); \
	fflush(stdout);\
  }
#endif
#define COUNT(name,test) {                                              \
    size_t notok = test;                                                \
    printf("%30s ... %s (%f%%)\n",name,notok ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)notok/(double)testN); \
	fflush(stdout);\
  }

  const size_t numSceneFlags = 64;

  RTCSceneFlags getSceneFlag(size_t i) 
  {
    int flag = 0;                               
    if (i & 1) flag |= RTC_SCENE_DYNAMIC;
    if (i & 2) flag |= RTC_SCENE_COMPACT;
    if (i & 4) flag |= RTC_SCENE_COHERENT;
    if (i & 8) flag |= RTC_SCENE_INCOHERENT;
    if (i & 16) flag |= RTC_SCENE_HIGH_QUALITY;
    if (i & 32) flag |= RTC_SCENE_ROBUST;
    return (RTCSceneFlags) flag;
  }

  RTCRay makeRay(const Vec3fa& org, const Vec3fa& dir) 
  {
    RTCRay ray;
    ray.org[0] = org.x; ray.org[1] = org.y; ray.org[2] = org.z;
    ray.dir[0] = dir.x; ray.dir[1] = dir.y; ray.dir[2] = dir.z;
    ray.tnear = 0.0f; ray.tfar = inf;
    ray.time = 0; ray.mask = -1;
    ray.geomID = ray.primID = ray.instID = -1;
    return ray;
  }

  RTCRay makeRay(const Vec3fa& org, const Vec3fa& dir, float tnear, float tfar) 
  {
    RTCRay ray;
    ray.org[0] = org.x; ray.org[1] = org.y; ray.org[2] = org.z;
    ray.dir[0] = dir.x; ray.dir[1] = dir.y; ray.dir[2] = dir.z;
    ray.tnear = tnear; ray.tfar = tfar;
    ray.time = 0; ray.mask = -1;
    ray.geomID = ray.primID = ray.instID = -1;
    return ray;
  }
  
  void setRay(RTCRay4& ray_o, int i, const RTCRay& ray_i)
  {
    ray_o.orgx[i] = ray_i.org[0];
    ray_o.orgy[i] = ray_i.org[1];
    ray_o.orgz[i] = ray_i.org[2];
    ray_o.dirx[i] = ray_i.dir[0];
    ray_o.diry[i] = ray_i.dir[1];
    ray_o.dirz[i] = ray_i.dir[2];
    ray_o.tnear[i] = ray_i.tnear;
    ray_o.tfar[i] = ray_i.tfar;
    ray_o.Ngx[i] = ray_i.Ng[0];
    ray_o.Ngy[i] = ray_i.Ng[1];
    ray_o.Ngz[i] = ray_i.Ng[2];
    ray_o.time[i] = ray_i.time;
    ray_o.mask[i] = ray_i.mask;
    ray_o.geomID[i] = ray_i.geomID;
    ray_o.primID[i] = ray_i.primID;
    ray_o.instID[i] = ray_i.instID;
  }

  void setRay(RTCRay8& ray_o, int i, const RTCRay& ray_i)
  {
    ray_o.orgx[i] = ray_i.org[0];
    ray_o.orgy[i] = ray_i.org[1];
    ray_o.orgz[i] = ray_i.org[2];
    ray_o.dirx[i] = ray_i.dir[0];
    ray_o.diry[i] = ray_i.dir[1];
    ray_o.dirz[i] = ray_i.dir[2];
    ray_o.tnear[i] = ray_i.tnear;
    ray_o.tfar[i] = ray_i.tfar;
    ray_o.Ngx[i] = ray_i.Ng[0];
    ray_o.Ngy[i] = ray_i.Ng[1];
    ray_o.Ngz[i] = ray_i.Ng[2];
    ray_o.time[i] = ray_i.time;
    ray_o.mask[i] = ray_i.mask;
    ray_o.geomID[i] = ray_i.geomID;
    ray_o.primID[i] = ray_i.primID;
    ray_o.instID[i] = ray_i.instID;
  }

  void setRay(RTCRay16& ray_o, int i, const RTCRay& ray_i)
  {
    ray_o.orgx[i] = ray_i.org[0];
    ray_o.orgy[i] = ray_i.org[1];
    ray_o.orgz[i] = ray_i.org[2];
    ray_o.dirx[i] = ray_i.dir[0];
    ray_o.diry[i] = ray_i.dir[1];
    ray_o.dirz[i] = ray_i.dir[2];
    ray_o.tnear[i] = ray_i.tnear;
    ray_o.tfar[i] = ray_i.tfar;
    ray_o.Ngx[i] = ray_i.Ng[0];
    ray_o.Ngy[i] = ray_i.Ng[1];
    ray_o.Ngz[i] = ray_i.Ng[2];
    ray_o.time[i] = ray_i.time;
    ray_o.mask[i] = ray_i.mask;
    ray_o.geomID[i] = ray_i.geomID;
    ray_o.primID[i] = ray_i.primID;
    ray_o.instID[i] = ray_i.instID;
  }

  RTCRay getRay(RTCRay4& ray_i, int i)
  {
    RTCRay ray_o;
    ray_o.org[0] = ray_i.orgx[i];
    ray_o.org[1] = ray_i.orgy[i];
    ray_o.org[2] = ray_i.orgz[i];
    ray_o.dir[0] = ray_i.dirx[i];
    ray_o.dir[1] = ray_i.diry[i];
    ray_o.dir[2] = ray_i.dirz[i];
    ray_o.tnear = ray_i.tnear[i];
    ray_o.tfar = ray_i.tfar[i];
    ray_o.Ng[0] = ray_i.Ngx[i];
    ray_o.Ng[1] = ray_i.Ngy[i];
    ray_o.Ng[2] = ray_i.Ngz[i];
    ray_o.time = ray_i.time[i];
    ray_o.mask = ray_i.mask[i];
    ray_o.geomID = ray_i.geomID[i];
    ray_o.primID = ray_i.primID[i];
    ray_o.instID = ray_i.instID[i];
    return ray_o;
  }

  RTCRay getRay(RTCRay8& ray_i, int i)
  {
    RTCRay ray_o;
    ray_o.org[0] = ray_i.orgx[i];
    ray_o.org[1] = ray_i.orgy[i];
    ray_o.org[2] = ray_i.orgz[i];
    ray_o.dir[0] = ray_i.dirx[i];
    ray_o.dir[1] = ray_i.diry[i];
    ray_o.dir[2] = ray_i.dirz[i];
    ray_o.tnear = ray_i.tnear[i];
    ray_o.tfar = ray_i.tfar[i];
    ray_o.Ng[0] = ray_i.Ngx[i];
    ray_o.Ng[1] = ray_i.Ngy[i];
    ray_o.Ng[2] = ray_i.Ngz[i];
    ray_o.time = ray_i.time[i];
    ray_o.mask = ray_i.mask[i];
    ray_o.geomID = ray_i.geomID[i];
    ray_o.primID = ray_i.primID[i];
    ray_o.instID = ray_i.instID[i];
    return ray_o;
  }

  RTCRay getRay(RTCRay16& ray_i, int i)
  {
    RTCRay ray_o;
    ray_o.org[0] = ray_i.orgx[i];
    ray_o.org[1] = ray_i.orgy[i];
    ray_o.org[2] = ray_i.orgz[i];
    ray_o.dir[0] = ray_i.dirx[i];
    ray_o.dir[1] = ray_i.diry[i];
    ray_o.dir[2] = ray_i.dirz[i];
    ray_o.tnear = ray_i.tnear[i];
    ray_o.tfar = ray_i.tfar[i];
    ray_o.Ng[0] = ray_i.Ngx[i];
    ray_o.Ng[1] = ray_i.Ngy[i];
    ray_o.Ng[2] = ray_i.Ngz[i];
    ray_o.time = ray_i.time[i];
    ray_o.mask = ray_i.mask[i];
    ray_o.geomID = ray_i.geomID[i];
    ray_o.primID = ray_i.primID[i];
    ray_o.instID = ray_i.instID[i];
    return ray_o;
  }

  void rtcIntersectN(RTCScene scene, RTCRay& ray, int N) 
  {
    switch (N) {
    case 1: {
      rtcIntersect(scene,ray); 
      break;
    }
#if !defined(__MIC__)
    case 4: {
      RTCRay4 ray4;
      for (size_t i=0; i<4; i++) setRay(ray4,i,ray);
      __aligned(16) int valid[4] = { -1,-1,-1,-1 };
      rtcIntersect4(valid,scene,ray4);
      ray = getRay(ray4,0);
      break;
    }
#endif
#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    case 8: {
      RTCRay8 ray8;
      for (size_t i=0; i<8; i++) setRay(ray8,i,ray);
      __aligned(32) int valid[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect8(valid,scene,ray8);
      ray = getRay(ray8,0);
      break;
    }
#endif
#if defined(__MIC__)
    case 16: {
      RTCRay16 ray16;
      for (size_t i=0; i<16; i++) setRay(ray16,i,ray);
      __aligned(64) int valid[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect16(valid,scene,ray16);
      ray = getRay(ray16,0);
      break;
    }
#endif
    default: break;
    }
  }

  void rtcOccludedN(RTCScene scene, RTCRay& ray, int N) 
  {
    switch (N) {
    case 1: {
      rtcOccluded(scene,ray); 
      break;
    }
#if !defined(__MIC__)
    case 4: {
      RTCRay4 ray4;
      for (size_t i=0; i<4; i++) setRay(ray4,i,ray);
      __aligned(16) int valid[4] = { -1,-1,-1,-1 };
      rtcOccluded4(valid,scene,ray4);
      ray.geomID = ray4.geomID[0];
      break;
    }
#endif
#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    case 8: {
      RTCRay8 ray8;
      for (size_t i=0; i<8; i++) setRay(ray8,i,ray);
      __aligned(32) int valid[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcOccluded8(valid,scene,ray8);
      ray.geomID = ray8.geomID[0];
      break;
    }
#endif
#if defined(__MIC__)
    case 16: {
      RTCRay16 ray16;
      for (size_t i=0; i<16; i++) setRay(ray16,i,ray);
      __aligned(64) int valid[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
      rtcOccluded16(valid,scene,ray16);
      ray.geomID = ray16.geomID[0];
      break;
    }
#endif
    default: break;
    }
  }

  static void parseCommandLine(int argc, char** argv)
  {
    for (int i=1; i<argc; i++)
    {
      std::string tag = argv[i];
      if (tag == "") return;

      /* rtcore configuration */
      else if (tag == "-rtcore" && i+1<argc) {
        g_rtcore = argv[++i];
      }

      /* skip unknown command line parameter */
      else {
        std::cerr << "unknown command line parameter: " << tag << " ";
        std::cerr << std::endl;
      }
    }
  }

  unsigned addPlane (RTCScene scene, RTCGeometryFlags flag, size_t num, const Vec3fa& p0, const Vec3fa& dx, const Vec3fa& dy)
  {
    unsigned mesh = rtcNewTriangleMesh (scene, flag, 2*num*num, (num+1)*(num+1));
    Vertex*   vertices  = (Vertex*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    for (size_t y=0; y<=num; y++) {
      for (size_t x=0; x<=num; x++) {
        Vec3fa p = p0+float(x)/float(num)*dx+float(y)/float(num)*dy;
        size_t i = y*(num+1)+x;
        vertices[i].x = p.x;
        vertices[i].y = p.y;
        vertices[i].z = p.z;
      }
    }
    for (size_t y=0; y<num; y++) {
      for (size_t x=0; x<num; x++) {
        size_t i = 2*y*num+2*x;
        size_t p00 = (y+0)*(num+1)+(x+0);
        size_t p01 = (y+0)*(num+1)+(x+1);
        size_t p10 = (y+1)*(num+1)+(x+0);
        size_t p11 = (y+1)*(num+1)+(x+1);
        triangles[i+0].v0 = p01; triangles[i+0].v1 = p00; triangles[i+0].v2 = p11;
        triangles[i+1].v0 = p10; triangles[i+1].v1 = p11; triangles[i+1].v2 = p00;
      }
    }
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    return mesh;
  }

  unsigned addSphere (RTCScene scene, RTCGeometryFlags flag, const Vec3fa& pos, const float r, size_t numPhi, size_t maxTriangles = -1, float motion = 0.0f)
  {
    /* create a triangulated sphere */
    size_t numTheta = 2*numPhi;
    size_t numTriangles = min(maxTriangles,2*numTheta*(numPhi-1));
    size_t numTimeSteps = motion == 0.0f ? 1 : 2;

    unsigned mesh = rtcNewTriangleMesh (scene, flag, numTriangles, numTheta*(numPhi+1),numTimeSteps);
    
    /* map triangle and vertex buffer */
    Vertex* vertices0 = NULL;
    Vertex* vertices1 = NULL;
    vertices0 = (Vertex*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
    if (numTimeSteps == 2) vertices1 = (Vertex*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);

    /* create sphere geometry */
    size_t tri = 0;
    const float rcpNumTheta = 1.0f/float(numTheta);
    const float rcpNumPhi   = 1.0f/float(numPhi);
    for (size_t phi=0; phi<=numPhi; phi++)
    {
      for (size_t theta=0; theta<numTheta; theta++)
      {
        const float phif   = phi*float(pi)*rcpNumPhi;
        const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
        Vertex* v = &vertices0[phi*numTheta+theta];
        const float cosThetaf = cos(thetaf);
        v->x = pos.x + r*sin(phif)*sin(thetaf);
        v->y = pos.y + r*cos(phif);
        v->z = pos.z + r*sin(phif)*cosThetaf;

        if (vertices1) {
          Vertex* v1 = &vertices1[phi*numTheta+theta];
          const float cosThetaf = cos(thetaf);
          v1->x = motion + pos.x + r*sin(phif)*sin(thetaf);
          v1->y = motion + pos.y + r*cos(phif);
          v1->z = motion + pos.z + r*sin(phif)*cosThetaf;
        }
      }
      if (phi == 0) continue;

      for (size_t theta=1; theta<=numTheta; theta++) 
      {
        int p00 = (phi-1)*numTheta+theta-1;
        int p01 = (phi-1)*numTheta+theta%numTheta;
        int p10 = phi*numTheta+theta-1;
        int p11 = phi*numTheta+theta%numTheta;
        
        if (phi > 1) {
          if (tri < numTriangles) {
            triangles[tri].v0 = p10; 
            triangles[tri].v1 = p00; 
            triangles[tri].v2 = p01; 
            tri++;
          }
        }
        
        if (phi < numPhi) {
          if (tri < numTriangles) {
            triangles[tri].v0 = p11; 
            triangles[tri].v1 = p10;
            triangles[tri].v2 = p01; 
            tri++;
          }
        }
      }
    }

    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
    if (numTimeSteps == 2) rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    return mesh;
  }

  struct Sphere
  {
    ALIGNED_CLASS;
  public:
    Vec3fa pos;
    float r;
  };

  void BoundsFunc(Sphere* sphere, size_t index, BBox3f* bounds_o)
  {
    bounds_o->lower.x = sphere->pos.x-sphere->r;
    bounds_o->lower.y = sphere->pos.y-sphere->r;
    bounds_o->lower.z = sphere->pos.z-sphere->r;
    bounds_o->upper.x = sphere->pos.x+sphere->r;
    bounds_o->upper.y = sphere->pos.y+sphere->r;
    bounds_o->upper.z = sphere->pos.z+sphere->r;
  }

  void IntersectFunc(void* ptr, RTCRay& ray, size_t item) {
  }

  void IntersectFunc4(const void* valid, void* ptr, RTCRay4& ray, size_t item) {
  }

  void IntersectFunc8(const void* valid, void* ptr, RTCRay8& ray, size_t item) {
  }

  void IntersectFunc16(const void* valid, void* ptr, RTCRay16& ray, size_t item) {
  }

  void OccludedFunc (void* ptr, RTCRay& ray, size_t item) {
  }

  void OccludedFunc4 (const void* valid, void* ptr, RTCRay4& ray, size_t item) {
  }

  void OccludedFunc8 (const void* valid, void* ptr, RTCRay8& ray, size_t item) {
  }

  void OccludedFunc16 (const void* valid, void* ptr, RTCRay16& ray, size_t item) {
  }

  unsigned addUserGeometryEmpty (RTCScene scene, const Vec3fa& pos, const float r)
  {
    BBox3f bounds(pos-Vec3fa(r),pos+Vec3fa(r));
    unsigned geom = rtcNewUserGeometry (scene,1);
    Sphere* sphere = new Sphere; // FIXME: get never deleted
    sphere->pos = pos; sphere->r = r;
    rtcSetBoundsFunction(scene,geom,(RTCBoundsFunc)BoundsFunc);
    rtcSetUserData(scene,geom,sphere);
    rtcSetIntersectFunction(scene,geom,IntersectFunc);
    rtcSetIntersectFunction4(scene,geom,IntersectFunc4);
    rtcSetIntersectFunction8(scene,geom,IntersectFunc8);
    rtcSetIntersectFunction16(scene,geom,&IntersectFunc16);
    rtcSetOccludedFunction(scene,geom,OccludedFunc);
    rtcSetOccludedFunction4(scene,geom,OccludedFunc4);
    rtcSetOccludedFunction8(scene,geom,OccludedFunc8);
    rtcSetOccludedFunction16(scene,geom,&OccludedFunc16);
    return geom;
  }

  BarrierSys g_barrier;
  volatile atomic_t g_atomic0;
  volatile atomic_t g_atomic1;

  void test_barrier_sys_thread(void* ptr) 
  {
    for (size_t i=0; i<10000; i++) 
    {
      atomic_add(&g_atomic0,+1);
      g_barrier.wait();
      atomic_add(&g_atomic1,+1);
      g_barrier.wait();
      atomic_add(&g_atomic0,-1);
      g_barrier.wait();
      atomic_add(&g_atomic1,-1);
      g_barrier.wait();
    }
  }

  bool test_barrier_sys ()
  {
    size_t numThreads = getNumberOfLogicalThreads();
#if defined (__MIC__)
    numThreads -= 4;
#endif
    g_barrier.init(numThreads);
    g_atomic0 = 0;
    g_atomic1 = 0;
    for (size_t i=1; i<numThreads; i++)
      g_threads.push_back(createThread(test_barrier_sys_thread,NULL,1000000,i));
    setAffinity(0);
    
    bool ok = true;
    for (size_t i=0; i<10000; i++) 
    {
      atomic_add(&g_atomic0,+1);
      g_barrier.wait();
      if (g_atomic0 != numThreads) ok = false;

      atomic_add(&g_atomic1,+1);
      g_barrier.wait();
      if (g_atomic1 != numThreads) ok = false;

      atomic_add(&g_atomic0,-1);
      g_barrier.wait();
      if (g_atomic0 != 0) ok = false;

      atomic_add(&g_atomic1,-1);
      g_barrier.wait();
      if (g_atomic1 != 0) ok = false;
    }

    for (size_t i=0; i<g_threads.size(); i++)
      join(g_threads[i]);

    g_threads.clear();
    return ok;
  }

  struct BarrierUsingCondition
  {
    __forceinline BarrierUsingCondition () 
      : count(0), barrierSize(0) {}
    
    __forceinline void init(size_t N) 
    {
      count = 0;
      barrierSize = N;
    }

    __forceinline void wait(int threadIndex)
    {
      mutex.lock();
      count++;

      if (count == barrierSize) {
        count = 0;
        cond.broadcast();
        mutex.unlock();
        return;
      }
     
      cond.wait(mutex);
      mutex.unlock();
      return;
    }

  public:
    MutexSys mutex;
    ConditionSys cond;
    volatile atomic_t count;
    volatile size_t barrierSize;
  };

  BarrierUsingCondition g_cond_barrier;

  void test_condition_sys_thread(void* ptr) 
  {
    for (size_t i=0; i<1000; i++) 
    {
      atomic_add(&g_atomic0,+1);
      g_cond_barrier.wait(1);
      atomic_add(&g_atomic1,+1);
      g_cond_barrier.wait(1);
      atomic_add(&g_atomic0,-1);
      g_cond_barrier.wait(1);
      atomic_add(&g_atomic1,-1);
      g_cond_barrier.wait(1);
    }
  }

  bool test_condition_sys ()
  {
    size_t numThreads = getNumberOfLogicalThreads();
#if defined (__MIC__)
    numThreads -= 4;
#endif
    g_cond_barrier.init(numThreads);
    g_atomic0 = 0;
    g_atomic1 = 0;
    for (size_t i=1; i<numThreads; i++)
      g_threads.push_back(createThread(test_condition_sys_thread,NULL,1000000,i));
    setAffinity(0);
    
    bool ok = true;
    for (size_t i=0; i<1000; i++) 
    {
      atomic_add(&g_atomic0,+1);
      g_cond_barrier.wait(0);
      if (g_atomic0 != numThreads) ok = false;

      atomic_add(&g_atomic1,+1);
      g_cond_barrier.wait(0);
      if (g_atomic1 != numThreads) ok = false;

      atomic_add(&g_atomic0,-1);
      g_cond_barrier.wait(0);
      if (g_atomic0 != 0) ok = false;

      atomic_add(&g_atomic1,-1);
      g_cond_barrier.wait(0);
      if (g_atomic1 != 0) ok = false;
    }

    for (size_t i=0; i<g_threads.size(); i++)
      join(g_threads[i]);

    g_threads.clear();
    return ok;
  }

  MutexSys g_mutex;
  size_t g_counter;

  void test_mutex_sys_thread(void* ptr) 
  {
    for (size_t i=0; i<10000; i++) 
    {
      g_mutex.lock();
      g_counter++;
      g_mutex.unlock();
    }
  }

  bool test_mutex_sys ()
  {
    size_t numThreads = getNumberOfLogicalThreads();
#if defined (__MIC__)
    numThreads -= 4;
#endif
    g_barrier.init(numThreads);
    g_counter = 0;
    for (size_t i=1; i<numThreads; i++) 
      g_threads.push_back(createThread(test_mutex_sys_thread,NULL,1000000,i));

    setAffinity(0);
    
    for (size_t i=0; i<10000; i++) 
    {
      g_mutex.lock();
      g_counter++;
      g_mutex.unlock();
    }

    for (size_t i=0; i<g_threads.size(); i++)
      join(g_threads[i]);

    g_threads.clear();
    return g_counter == 10000*numThreads;
  }

  bool rtcore_empty(RTCSceneFlags flags)
  {
    RTCScene scene = rtcNewScene(flags,aflags);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  bool rtcore_dynamic_flag(RTCSceneFlags sceneFlag, RTCGeometryFlags geomFlag)
  {
    RTCScene scene = rtcNewScene(sceneFlag,aflags);
    AssertNoError();
    unsigned geom = rtcNewTriangleMesh (scene, geomFlag, 0, 0);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  bool rtcore_static_scene()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    unsigned geom = addSphere(scene,RTC_GEOMETRY_STATIC,zero,1.0f,50);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
#if !defined(__EXIT_ON_ERROR__)
    rtcCommit (scene); // cannot commit static scene twice
    AssertAnyError();
    rtcDisable(scene,geom); // static scene cannot get modified anymore after commit
    AssertAnyError();
#endif
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  bool rtcore_deformable_geometry()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    unsigned geom = addSphere(scene,RTC_GEOMETRY_DEFORMABLE,zero,1.0f,50);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
#if !defined(__EXIT_ON_ERROR__)
    rtcMapBuffer(scene,geom,RTC_INDEX_BUFFER);
    AssertError(RTC_INVALID_OPERATION); // cannot modify index buffer of deformable geometry anymore after commit
#endif
    rtcMapBuffer(scene,geom,RTC_VERTEX_BUFFER);
    AssertNoError();
#if !defined(__EXIT_ON_ERROR__)
    rtcUnmapBuffer(scene,geom,RTC_INDEX_BUFFER);
    AssertError(RTC_INVALID_OPERATION); // cannot modify index buffer of deformable geometry anymore after commit
#endif
    rtcUnmapBuffer(scene,geom,RTC_VERTEX_BUFFER);
    AssertNoError();
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  bool rtcore_unmapped_before_commit()
  {
#if !defined(__EXIT_ON_ERROR__)
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    unsigned geom0 = addSphere(scene,RTC_GEOMETRY_STATIC,zero,1.0f,50);
    unsigned geom1 = addSphere(scene,RTC_GEOMETRY_STATIC,zero,1.0f,50);
    AssertNoError();
    rtcMapBuffer(scene,geom0,RTC_INDEX_BUFFER);
    rtcMapBuffer(scene,geom0,RTC_VERTEX_BUFFER);
    AssertNoError();
    rtcCommit (scene);
    AssertError(RTC_INVALID_OPERATION); // error, buffers still mapped
    rtcDeleteScene (scene);
    AssertNoError();
#endif
    return true;
  }

  bool rtcore_buffer_stride()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    unsigned geom = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, 16, 16);
    AssertNoError();
    char* indexBuffer  = (char*) alignedMalloc(8+16*6*sizeof(int));
    char* vertexBuffer = (char*) alignedMalloc(12+16*9*sizeof(float)+4);


#if !defined(__MIC__)    
    rtcSetBuffer(scene,geom,RTC_INDEX_BUFFER,indexBuffer,1,3*sizeof(int));
    AssertError(RTC_INVALID_OPERATION);
    rtcSetBuffer(scene,geom,RTC_VERTEX_BUFFER,vertexBuffer,1,3*sizeof(float));
    AssertError(RTC_INVALID_OPERATION);

    rtcSetBuffer(scene,geom,RTC_INDEX_BUFFER,indexBuffer,0,3*sizeof(int)+3);
    AssertError(RTC_INVALID_OPERATION);
    rtcSetBuffer(scene,geom,RTC_VERTEX_BUFFER,vertexBuffer,0,3*sizeof(float)+3);
    AssertError(RTC_INVALID_OPERATION);

    rtcSetBuffer(scene,geom,RTC_INDEX_BUFFER,indexBuffer,0,3*sizeof(int));
    AssertNoError();
    rtcSetBuffer(scene,geom,RTC_VERTEX_BUFFER,vertexBuffer,0,3*sizeof(float));
    AssertNoError();

    rtcSetBuffer(scene,geom,RTC_INDEX_BUFFER,indexBuffer,8,6*sizeof(int));
    AssertNoError();
    rtcSetBuffer(scene,geom,RTC_VERTEX_BUFFER,vertexBuffer,12,9*sizeof(float));
    AssertNoError();


    rtcSetBuffer(scene,geom,RTC_INDEX_BUFFER,indexBuffer,0,3*sizeof(int));
    AssertNoError();
#endif

    rtcSetBuffer(scene,geom,RTC_VERTEX_BUFFER,vertexBuffer,0,4*sizeof(float));
    AssertNoError();

    rtcDeleteScene (scene);
    AssertNoError();

    alignedFree(indexBuffer);
    alignedFree(vertexBuffer);
    return true;
  }

  bool rtcore_dynamic_enable_disable()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    unsigned geom0 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(-1,0,-1),1.0f,50);
    unsigned geom1 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(-1,0,+1),1.0f,50);
    unsigned geom2 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(+1,0,-1),1.0f,50);
    unsigned geom3 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(+1,0,+1),1.0f,50);
    AssertNoError();

    for (size_t i=0; i<16; i++) 
    {
      bool enabled0 = i & 1, enabled1 = i & 2, enabled2 = i & 4, enabled3 = i & 8;
      if (enabled0) rtcEnable(scene,geom0); else rtcDisable(scene,geom0); AssertNoError();
      if (enabled1) rtcEnable(scene,geom1); else rtcDisable(scene,geom1); AssertNoError();
      if (enabled2) rtcEnable(scene,geom2); else rtcDisable(scene,geom2); AssertNoError();
      if (enabled3) rtcEnable(scene,geom3); else rtcDisable(scene,geom3); AssertNoError();
      rtcCommit (scene);
      AssertNoError();
      {
        RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0));
        RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); 
        RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); 
        RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); 
        rtcIntersect(scene,ray0);
        rtcIntersect(scene,ray1);
        rtcIntersect(scene,ray2);
        rtcIntersect(scene,ray3);
        bool ok0 = enabled0 ? ray0.geomID == 0 : ray0.geomID == -1;
        bool ok1 = enabled1 ? ray1.geomID == 1 : ray1.geomID == -1;
        bool ok2 = enabled2 ? ray2.geomID == 2 : ray2.geomID == -1;
        bool ok3 = enabled3 ? ray3.geomID == 3 : ray3.geomID == -1;
        if (!ok0 || !ok1 || !ok2 || !ok3) return false;
      }
    }
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  void move_mesh(RTCScene scene, unsigned mesh, size_t numVertices, Vec3fa& pos) 
  {
    Vec3fa* vertices = (Vec3fa*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    for (size_t i=0; i<numVertices; i++) vertices[i] += pos;
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER);
    rtcUpdate(scene,mesh);
  }
  
  bool rtcore_update(RTCGeometryFlags flags)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    size_t numPhi = 50;
    size_t numVertices = 2*numPhi*(numPhi+1);
    unsigned geom0 = addSphere(scene,flags,Vec3fa(-1,0,-1),1.0f,numPhi);
    unsigned geom1 = addSphere(scene,flags,Vec3fa(-1,0,+1),1.0f,numPhi);
    unsigned geom2 = addSphere(scene,flags,Vec3fa(+1,0,-1),1.0f,numPhi);
    unsigned geom3 = addSphere(scene,flags,Vec3fa(+1,0,+1),1.0f,numPhi);
    Vec3fa pos0 = zero;
    Vec3fa pos1 = zero;
    Vec3fa pos2 = zero;
    Vec3fa pos3 = zero;
    AssertNoError();
    
    for (size_t i=0; i<16; i++) 
    {
      bool move0 = i & 1, move1 = i & 2, move2 = i & 4, move3 = i & 8;
      Vec3fa ds(20,0,20);
      if (move0) { move_mesh(scene,geom0,numVertices,ds); pos0 += ds; }
      if (move1) { move_mesh(scene,geom1,numVertices,ds); pos1 += ds; }
      if (move2) { move_mesh(scene,geom2,numVertices,ds); pos2 += ds; }
      if (move3) { move_mesh(scene,geom3,numVertices,ds); pos3 += ds; }
      rtcCommit (scene);
      AssertNoError();
      {
        RTCRay ray0 = makeRay(pos0+Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); 
        RTCRay ray1 = makeRay(pos1+Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); 
        RTCRay ray2 = makeRay(pos2+Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); 
        RTCRay ray3 = makeRay(pos3+Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); 
        rtcIntersect(scene,ray0);
        rtcIntersect(scene,ray1);
        rtcIntersect(scene,ray2);
        rtcIntersect(scene,ray3);
        if (ray0.geomID != 0 || 
            ray1.geomID != 1 || 
            ray2.geomID != 2 || 
            ray3.geomID != 3) return false;

#if !defined(__MIC__)
        RTCRay4 ray4; 
        setRay(ray4,0,ray0);
        setRay(ray4,1,ray1);
        setRay(ray4,2,ray2);
        setRay(ray4,3,ray3);
        __aligned(16) int valid4[4] = { -1,-1,-1,-1 };
        rtcIntersect4(valid4,scene,ray4);
        if (ray4.geomID[0] != 0 || 
            ray4.geomID[1] != 1 || 
            ray4.geomID[2] != 2 || 
            ray4.geomID[3] != 3) return false;
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
        if (has_feature(AVX)) 
        {
          RTCRay8 ray8; 
          setRay(ray8,0,ray0);
          setRay(ray8,1,ray1);
          setRay(ray8,2,ray2);
          setRay(ray8,3,ray3);
          __aligned(32) int valid8[8] = { -1,-1,-1,-1, 0, 0, 0, 0 };
          rtcIntersect8(valid8,scene,ray8);
          if (ray8.geomID[0] != 0 || 
              ray8.geomID[1] != 1 || 
              ray8.geomID[2] != 2 || 
              ray8.geomID[3] != 3) return false;
        }
#endif

#if defined(__MIC__)
        RTCRay16 ray16; 
        setRay(ray16,0,ray0);
        setRay(ray16,1,ray1);
        setRay(ray16,2,ray2);
        setRay(ray16,3,ray3);
        __aligned(64) int valid16[16] = { -1,-1,-1,-1,+0,+0,+0,+0, 
                            +0,+0,+0,+0,+0,+0,+0,+0 };
        rtcIntersect16(valid16,scene,ray16);
        if (ray16.geomID[0] != 0 || 
            ray16.geomID[1] != 1 || 
            ray16.geomID[2] != 2 || 
            ray16.geomID[3] != 3) return false;
#endif
      }
    }
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  bool rtcore_ray_masks_intersect(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;

    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom0 = addSphere(scene,gflags,Vec3fa(-1,0,-1),1.0f,50);
    unsigned geom1 = addSphere(scene,gflags,Vec3fa(-1,0,+1),1.0f,50);
    unsigned geom2 = addSphere(scene,gflags,Vec3fa(+1,0,-1),1.0f,50);
    unsigned geom3 = addSphere(scene,gflags,Vec3fa(+1,0,+1),1.0f,50);
    rtcSetMask(scene,geom0,1);
    rtcSetMask(scene,geom1,2);
    rtcSetMask(scene,geom2,4);
    rtcSetMask(scene,geom3,8);
    rtcCommit (scene);

    for (size_t i=0; i<16; i++) 
    {
      int mask0 = i;
      int mask1 = i+1;
      int mask2 = i+2;
      int mask3 = i+3;

      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;
	rtcIntersect(scene,ray0);
	rtcIntersect(scene,ray1);
	rtcIntersect(scene,ray2);
	rtcIntersect(scene,ray3);
	bool ok0 = mask0 & 1 ? ray0.geomID == 0 : ray0.geomID == -1;
	bool ok1 = mask1 & 2 ? ray1.geomID == 1 : ray1.geomID == -1;
	bool ok2 = mask2 & 4 ? ray2.geomID == 2 : ray2.geomID == -1;
	bool ok3 = mask3 & 8 ? ray3.geomID == 3 : ray3.geomID == -1;
	if (!ok0 || !ok1 || !ok2 || !ok3) passed = false;
      }

#if !defined(__MIC__)
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay4 ray4;
	setRay(ray4,0,ray0);
	setRay(ray4,1,ray1);
	setRay(ray4,2,ray2);
	setRay(ray4,3,ray3);
	__aligned(16) int valid4[4] = { -1,-1,-1,-1 };
	rtcIntersect4(valid4,scene,ray4);
	bool ok4a = mask0 & 1 ? ray4.geomID[0] == 0 : ray4.geomID[0] == -1;
	bool ok4b = mask1 & 2 ? ray4.geomID[1] == 1 : ray4.geomID[1] == -1;
	bool ok4c = mask2 & 4 ? ray4.geomID[2] == 2 : ray4.geomID[2] == -1;
	bool ok4d = mask3 & 8 ? ray4.geomID[3] == 3 : ray4.geomID[3] == -1;
	if (!ok4a || !ok4b || !ok4c || !ok4d) passed = false;
      }

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
      if (has_feature(AVX))
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay8 ray8;
	setRay(ray8,0,ray0);
	setRay(ray8,1,ray1);
	setRay(ray8,2,ray2);
	setRay(ray8,3,ray3);
	__aligned(32) int valid8[8] = { -1,-1,-1,-1,0,0,0,0 };
	rtcIntersect8(valid8,scene,ray8);
	bool ok8a = mask0 & 1 ? ray8.geomID[0] == 0 : ray8.geomID[0] == -1;
	bool ok8b = mask1 & 2 ? ray8.geomID[1] == 1 : ray8.geomID[1] == -1;
	bool ok8c = mask2 & 4 ? ray8.geomID[2] == 2 : ray8.geomID[2] == -1;
	bool ok8d = mask3 & 8 ? ray8.geomID[3] == 3 : ray8.geomID[3] == -1;
	if (!ok8a || !ok8b || !ok8c || !ok8d) passed = false;
      }
#endif

#endif

#if defined(__MIC__)
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay16 ray16;
	setRay(ray16,0,ray0);
	setRay(ray16,1,ray1);
	setRay(ray16,2,ray2);
	setRay(ray16,3,ray3);
	__aligned(64) int valid16[16] = { -1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0 };
	rtcIntersect16(valid16,scene,ray16);
	bool ok16a = mask0 & 1 ? ray16.geomID[0] == 0 : ray16.geomID[0] == -1;
	bool ok16b = mask1 & 2 ? ray16.geomID[1] == 1 : ray16.geomID[1] == -1;
	bool ok16c = mask2 & 4 ? ray16.geomID[2] == 2 : ray16.geomID[2] == -1;
	bool ok16d = mask3 & 8 ? ray16.geomID[3] == 3 : ray16.geomID[3] == -1;
	if (!ok16a || !ok16b || !ok16c || !ok16d) passed = false;
      }
#endif


    }
    rtcDeleteScene (scene);
    return passed;
  }

  bool rtcore_ray_masks_occluded(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;
    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom0 = addSphere(scene,gflags,Vec3fa(-1,0,-1),1.0f,50);
    unsigned geom1 = addSphere(scene,gflags,Vec3fa(-1,0,+1),1.0f,50);
    unsigned geom2 = addSphere(scene,gflags,Vec3fa(+1,0,-1),1.0f,50);
    unsigned geom3 = addSphere(scene,gflags,Vec3fa(+1,0,+1),1.0f,50);
    rtcSetMask(scene,geom0,1);
    rtcSetMask(scene,geom1,2);
    rtcSetMask(scene,geom2,4);
    rtcSetMask(scene,geom3,8);
    rtcCommit (scene);

    for (size_t i=0; i<16; i++) 
    {
      int mask0 = i;
      int mask1 = i+1;
      int mask2 = i+2;
      int mask3 = i+3;

      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;
	rtcOccluded(scene,ray0);
	rtcOccluded(scene,ray1);
	rtcOccluded(scene,ray2);
	rtcOccluded(scene,ray3);
	bool ok0 = mask0 & 1 ? ray0.geomID == 0 : ray0.geomID == -1;
	bool ok1 = mask1 & 2 ? ray1.geomID == 0 : ray1.geomID == -1;
	bool ok2 = mask2 & 4 ? ray2.geomID == 0 : ray2.geomID == -1;
	bool ok3 = mask3 & 8 ? ray3.geomID == 0 : ray3.geomID == -1;

	if (!ok0 || !ok1 || !ok2 || !ok3) passed = false;
      }

#if !defined(__MIC__)
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay4 ray4;
	setRay(ray4,0,ray0);
	setRay(ray4,1,ray1);
	setRay(ray4,2,ray2);
	setRay(ray4,3,ray3);
	__aligned(16) int valid4[4] = { -1,-1,-1,-1 };
	rtcOccluded4(valid4,scene,ray4);
	bool ok4a = mask0 & 1 ? ray4.geomID[0] == 0 : ray4.geomID[0] == -1;
	bool ok4b = mask1 & 2 ? ray4.geomID[1] == 0 : ray4.geomID[1] == -1;
	bool ok4c = mask2 & 4 ? ray4.geomID[2] == 0 : ray4.geomID[2] == -1;
	bool ok4d = mask3 & 8 ? ray4.geomID[3] == 0 : ray4.geomID[3] == -1;
	if (!ok4a || !ok4b || !ok4c || !ok4d) passed = false;
      }

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
      if (has_feature(AVX)) 
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay8 ray8;
	setRay(ray8,0,ray0);
	setRay(ray8,1,ray1);
	setRay(ray8,2,ray2);
	setRay(ray8,3,ray3);
	__aligned(32) int valid8[8] = { -1,-1,-1,-1,0,0,0,0 };
	rtcOccluded8(valid8,scene,ray8);
	bool ok8a = mask0 & 1 ? ray8.geomID[0] == 0 : ray8.geomID[0] == -1;
	bool ok8b = mask1 & 2 ? ray8.geomID[1] == 0 : ray8.geomID[1] == -1;
	bool ok8c = mask2 & 4 ? ray8.geomID[2] == 0 : ray8.geomID[2] == -1;
	bool ok8d = mask3 & 8 ? ray8.geomID[3] == 0 : ray8.geomID[3] == -1;
	if (!ok8a || !ok8b || !ok8c || !ok8d) passed = false;
      }
#endif

#endif

#if defined(__MIC__)
      {
	RTCRay ray0 = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(Vec3fa(-1,10,+1),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(Vec3fa(+1,10,-1),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(Vec3fa(+1,10,+1),Vec3fa(0,-1,0)); ray3.mask = mask3;

	RTCRay16 ray16;
	setRay(ray16,0,ray0);
	setRay(ray16,1,ray1);
	setRay(ray16,2,ray2);
	setRay(ray16,3,ray3);
	__aligned(64) int valid16[16] = { -1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0 };

	rtcOccluded16(valid16,scene,ray16);

	bool ok16a = mask0 & 1 ? ray16.geomID[0] == 0 : ray16.geomID[0] == -1;
	bool ok16b = mask1 & 2 ? ray16.geomID[1] == 0 : ray16.geomID[1] == -1;
	bool ok16c = mask2 & 4 ? ray16.geomID[2] == 0 : ray16.geomID[2] == -1;
	bool ok16d = mask3 & 8 ? ray16.geomID[3] == 0 : ray16.geomID[3] == -1;
	if (!ok16a || !ok16b || !ok16c || !ok16d) passed = false;
      }

#endif
    }
    rtcDeleteScene (scene);
    return passed;
  }
  
  void rtcore_ray_masks_all()
  {
    printf("%30s ... ","ray_masks");
    bool passed = true;
    for (int i=0; i<numSceneFlags; i++) 
    {
      RTCSceneFlags flag = getSceneFlag(i);
      bool ok0 = rtcore_ray_masks_intersect(flag,RTC_GEOMETRY_STATIC);
      if (ok0) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok0;
      bool ok1 = rtcore_ray_masks_occluded(flag,RTC_GEOMETRY_STATIC);
      if (ok1) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok1;
    }
    printf(" %s\n",passed ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m");
	fflush(stdout);
  }

  void intersectionFilter1(void* ptr, RTCRay& ray) 
  {
    if ((size_t)ptr != 123) 
      return;

    if (ray.primID & 2) 
      ray.geomID = -1;
  }

  void intersectionFilter4(const void* valid_i, void* ptr, RTCRay4& ray) 
  {
    if ((size_t)ptr != 123) 
      return;

    int* valid = (int*)valid_i;
    for (size_t i=0; i<4; i++)
      if (valid[i] == -1)
        if (ray.primID[i] & 2) 
          ray.geomID[i] = -1;
  }

  void intersectionFilter8(const void* valid_i, void* ptr, RTCRay8& ray) 
  {
    if ((size_t)ptr != 123) 
      return;

    int* valid = (int*)valid_i;
    for (size_t i=0; i<8; i++)
      if (valid[i] == -1)
        if (ray.primID[i] & 2) 
          ray.geomID[i] = -1;
  }

  void intersectionFilter16(const void* valid_i, void* ptr, RTCRay16& ray) 
  {
    if ((size_t)ptr != 123) 
      return;

    unsigned int valid = *(unsigned int*)valid_i;
    for (size_t i=0; i<16; i++)
	if (valid & ((unsigned int)1 << i))
	  if (ray.primID[i] & 2) 
	    ray.geomID[i] = -1;
  }

  bool rtcore_filter_intersect(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;

    RTCScene scene = rtcNewScene(sflags,aflags);
    Vec3fa p0(-0.75f,-0.25f,-10.0f), dx(4,0,0), dy(0,4,0);
    int geom0 = addPlane (scene, gflags, 4, p0, dx, dy);
    rtcSetUserData(scene,geom0,(void*)123);
    rtcSetIntersectionFilterFunction(scene,geom0,intersectionFilter1);
    rtcSetIntersectionFilterFunction4(scene,geom0,intersectionFilter4);
    rtcSetIntersectionFilterFunction8(scene,geom0,intersectionFilter8);
    rtcSetIntersectionFilterFunction16(scene,geom0,intersectionFilter16);
    rtcCommit (scene);
    
    for (size_t iy=0; iy<4; iy++) 
    {
      for (size_t ix=0; ix<4; ix++) 
      {
        int primID = 2*(iy*4+ix);
        {
          RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));
          rtcIntersect(scene,ray0);
          bool ok0 = (primID & 2) ? (ray0.geomID == -1) : (ray0.geomID == 0);
          if (!ok0) passed = false;
        }

#if !defined(__MIC__)
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay4 ray4;
	setRay(ray4,0,ray0);
	__aligned(16) int valid4[4] = { -1,0,0,0 };
	rtcIntersect4(valid4,scene,ray4);
        bool ok0 = (primID & 2) ? (ray4.geomID[0] == -1) : (ray4.geomID[0] == 0);
        if (!ok0) passed = false;
      }

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
      if (has_feature(AVX))
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay8 ray8;
	setRay(ray8,0,ray0);
	__aligned(32) int valid8[8] = { -1,0,0,0,0,0,0,0 };
	rtcIntersect8(valid8,scene,ray8);
        bool ok0 = (primID & 2) ? (ray8.geomID[0] == -1) : (ray8.geomID[0] == 0);
        if (!ok0) passed = false;
      }
#endif

#endif

#if defined(__MIC__) && 1
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay16 ray16;
	setRay(ray16,0,ray0);
	__aligned(64) int valid16[16] = { -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	rtcIntersect16(valid16,scene,ray16);
        bool ok0 = (primID & 2) ? (ray16.geomID[0] == -1) : (ray16.geomID[0] == 0);
        if (!ok0) passed = false;
      }
#endif
      }
    }
    rtcDeleteScene (scene);
    return passed;
  }

  bool rtcore_filter_occluded(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;

    RTCScene scene = rtcNewScene(sflags,aflags);
    Vec3fa p0(-0.75f,-0.25f,-10.0f), dx(4,0,0), dy(0,4,0);
    int geom0 = addPlane (scene, gflags, 4, p0, dx, dy);
    rtcSetUserData(scene,geom0,(void*)123);
    rtcSetOcclusionFilterFunction(scene,geom0,intersectionFilter1);
    rtcSetOcclusionFilterFunction4(scene,geom0,intersectionFilter4);
    rtcSetOcclusionFilterFunction8(scene,geom0,intersectionFilter8);
    rtcSetOcclusionFilterFunction16(scene,geom0,intersectionFilter16);
    rtcCommit (scene);
    
    for (size_t iy=0; iy<4; iy++) 
    {
      for (size_t ix=0; ix<4; ix++) 
      {
        int primID = 2*(iy*4+ix);
        {
          RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));
          rtcOccluded(scene,ray0);
          bool ok0 = (primID & 2) ? (ray0.geomID == -1) : (ray0.geomID == 0);
          if (!ok0) passed = false;
        }

#if !defined(__MIC__)
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay4 ray4;
	setRay(ray4,0,ray0);
	__aligned(16) int valid4[4] = { -1,0,0,0 };
	rtcOccluded4(valid4,scene,ray4);
        bool ok0 = (primID & 2) ? (ray4.geomID[0] == -1) : (ray4.geomID[0] == 0);
        if (!ok0) passed = false;
      }

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
      if (has_feature(AVX))
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay8 ray8;
	setRay(ray8,0,ray0);
	__aligned(32) int valid8[8] = { -1,0,0,0,0,0,0,0 };
	rtcOccluded8(valid8,scene,ray8);
        bool ok0 = (primID & 2) ? (ray8.geomID[0] == -1) : (ray8.geomID[0] == 0);
        if (!ok0) passed = false;
      }
#endif

#endif

#if defined(__MIC__)
      {
        RTCRay ray0 = makeRay(Vec3fa(float(ix),float(iy),0.0f),Vec3fa(0,0,-1));

	RTCRay16 ray16;
	setRay(ray16,0,ray0);
	__aligned(64) int valid16[16] = { -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	rtcOccluded16(valid16,scene,ray16);
        bool ok0 = (primID & 2) ? (ray16.geomID[0] == -1) : (ray16.geomID[0] == 0);
        if (!ok0) passed = false;
      }
#endif
      }
    }
    rtcDeleteScene (scene);
    return passed;
  }

  void rtcore_filter_all()
  {
    printf("%30s ... ","intersection_filter");
    bool passed = true;
    for (int i=0; i<numSceneFlags; i++) 
    {
      RTCSceneFlags flag = getSceneFlag(i);
      bool ok0 = rtcore_filter_intersect(flag,RTC_GEOMETRY_STATIC);
      if (ok0) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok0;
      bool ok1 = rtcore_filter_occluded(flag,RTC_GEOMETRY_STATIC);
      if (ok1) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok1;
    }
    printf(" %s\n",passed ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m");
    fflush(stdout);
  }

  bool rtcore_packet_write_test(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;
    RTCScene scene = rtcNewScene(sflags,aflags);

    addSphere(scene,gflags,Vec3fa(-1,0,-1),1.0f,50,-1,0.0f);

#if !defined(__MIC__)
    addSphere(scene,gflags,Vec3fa(+1,0,+1),1.0f,50,-1,0.1f);
#endif
    rtcCommit (scene);
    
    for (size_t i=0; i<4; i++) 
    {
      RTCRay ray = makeRay(Vec3fa(-1,10,-1),Vec3fa(0,-1,0));

#if !defined(__MIC__)
      RTCRay4 ray4; 
      memset(&ray4,-1,sizeof(RTCRay4));
      setRay(ray4,i,ray);
      __aligned(16) int valid4[4] = { 0,0,0,0 };
      valid4[i] = -1;
      rtcOccluded4(valid4,scene,ray4);
      rtcIntersect4(valid4,scene,ray4);
      
      for (int j=0; j<sizeof(RTCRay4)/4; j++) {
        if ((j%4) == i) continue;
        passed &= ((int*)&ray4)[j] == -1;
      }
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
      if (has_feature(AVX)) {
        RTCRay8 ray8; 
        memset(&ray8,-1,sizeof(RTCRay8));
        setRay(ray8,i,ray);
        __aligned(32) int valid8[8] = { 0,0,0,0,0,0,0,0 };
        valid8[i] = -1;
        rtcOccluded8(valid8,scene,ray8);
        rtcIntersect8(valid8,scene,ray8);
        
        for (int j=0; j<sizeof(RTCRay8)/4; j++) {
          if ((j%8) == i) continue;
          passed &= ((int*)&ray8)[j] == -1;
        }
      }
#endif

#if defined(__MIC__)
      __aligned(64) RTCRay16 ray16; 
      memset(&ray16,-1,sizeof(RTCRay16));
      setRay(ray16,i,ray);
      __aligned(64) int valid16[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
      valid16[i] = -1;
      rtcOccluded16(valid16,scene,ray16);
      rtcIntersect16(valid16,scene,ray16);
      
      for (int j=0; j<sizeof(RTCRay16)/4; j++) {
        if ((j%16) == i) continue;
        passed &= ((int*)&ray16)[j] == -1;
      }
#endif
    }
    return passed;
  }

  void rtcore_packet_write_test_all()
  {
    printf("%30s ... ","packet_write_test");
    bool passed = true;
    for (int i=0; i<numSceneFlags; i++) 
    {
      RTCSceneFlags flag = getSceneFlag(i);
      bool ok0 = rtcore_packet_write_test(flag,RTC_GEOMETRY_STATIC);
      if (ok0) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok0;
    }
    printf(" %s\n",passed ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m");
	  fflush(stdout);
  }

  void rtcore_watertight_sphere1(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(pos,0.0f,0.0f),2.0f,1000);
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(Vec3fa(pos,0.0f,0.0f)+org,dir); 
      rtcIntersect(scene,ray);
      numFailures += ray.primID == -1;
    }
    rtcDeleteScene (scene);

    printf("%30s ... %s (%f%%)\n","watertight_sphere1",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }
  
  void rtcore_watertight_sphere4(float pos)
  {
    RTCScene scene = rtcNewScene(RTCSceneFlags(RTC_SCENE_STATIC | RTC_SCENE_ROBUST),aflags);
    unsigned geom = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(pos,0.0f,0.0f),2.0f,1000);
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=4) {
      RTCRay4 ray4;
      for (size_t j=0; j<4; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos,0.0f,0.0f)+org,dir); 
        setRay(ray4,j,ray);
      }
      __aligned(16) int valid[4] = { -1,-1,-1,-1 };
      rtcIntersect4(valid,scene,ray4);
      for (size_t j=0; j<4; j++)
        numFailures += ray4.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_sphere4",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }

  void rtcore_watertight_sphere8(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(pos,0.0f,0.0f),2.0f,1000);
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=8) {
      RTCRay8 ray8;
      for (size_t j=0; j<8; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos,0.0f,0.0f)+org,dir); 
        setRay(ray8,j,ray);
      }
      __aligned(32) int valid[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect8(valid,scene,ray8);
      for (size_t j=0; j<8; j++)
        numFailures += ray8.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_sphere8",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }

  void rtcore_watertight_sphere16(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(pos,0.0f,0.0f),2.0f,1000);
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=16) {
      RTCRay16 ray16;
      for (size_t j=0; j<16; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos,0.0f,0.0f)+org,dir); 
        setRay(ray16,j,ray);
      }
      __aligned(64) int valid[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect16(valid,scene,ray16);
      for (size_t j=0; j<16; j++)
        numFailures += ray16.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_sphere16",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }
  
  void rtcore_watertight_plane1(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,1000,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i++) {
      Vec3fa org(drand48()-0.5f,drand48()-0.5f,drand48()-0.5f);
      Vec3fa dir(1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(Vec3fa(pos-3.0f,0.0f,0.0f),dir); 
      rtcIntersect(scene,ray);
      numFailures += ray.primID == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_plane1",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	fflush(stdout);
  }

  void rtcore_watertight_plane4(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,1000,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=4) {
      RTCRay4 ray4;
      for (size_t j=0; j<4; j++) {
        Vec3fa org(drand48()-0.5f,drand48()-0.5f,drand48()-0.5f);
        Vec3fa dir(1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos-3.0f,0.0f,0.0f),dir); 
        setRay(ray4,j,ray);
      }
      __aligned(16) int valid[4] = { -1,-1,-1,-1 };
      rtcIntersect4(valid,scene,ray4);
      for (size_t j=0; j<4; j++)
        numFailures += ray4.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_plane4",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
  	fflush(stdout);
  }

  void rtcore_watertight_plane8(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,1000,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=8) {
      RTCRay8 ray8;
      for (size_t j=0; j<8; j++) {
        Vec3fa org(drand48()-0.5f,drand48()-0.5f,drand48()-0.5f);
        Vec3fa dir(1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos-3.0f,0.0f,0.0f),dir); 
        setRay(ray8,j,ray);
      }
      __aligned(32) int valid[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect8(valid,scene,ray8);
      for (size_t j=0; j<8; j++)
        numFailures += ray8.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_plane8",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }

  void rtcore_watertight_plane16(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,1000,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
    rtcCommit (scene);
    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=16) {
      RTCRay16 ray16;
      for (size_t j=0; j<16; j++) {
        Vec3fa org(drand48()-0.5f,drand48()-0.5f,drand48()-0.5f);
        Vec3fa dir(1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(Vec3fa(pos-3.0f,0.0f,0.0f),dir); 
        setRay(ray16,j,ray);
      }
      __aligned(64) int valid[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect16(valid,scene,ray16);
      for (size_t j=0; j<16; j++)
        numFailures += ray16.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    printf("%30s ... %s (%f%%)\n","watertight_plane16",
           numFailures ? "\033[31m[FAILED]\033[0m" : "\033[32m[PASSED]\033[0m", 100.0f*(double)numFailures/(double)testN);
	  fflush(stdout);
  }

  void rtcore_nan(const char* name, RTCSceneFlags sflags, RTCGeometryFlags gflags, int N)
  {
    size_t count = 1000/N;
    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom = addSphere(scene,gflags,zero,2.0f,100);
    rtcCommit (scene);
    size_t numFailures = 0;
    //size_t c0 = __rdtsc();
    double c0 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c1 = __rdtsc();
    double c1 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(nan),dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c2 = __rdtsc();
    double c2 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(nan),dir+Vec3fa(nan)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c3 = __rdtsc();
    double c3 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir,nan,nan); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c4 = __rdtsc();
    double c4 = getSeconds();

    double d1 = c1-c0;
    double d2 = c2-c1;
    double d3 = c3-c2;
    double d4 = c4-c3;
    rtcDeleteScene (scene);

    bool ok = (d2 < 2.5*d1) && (d3 < 2.5*d1) && (d4 < 2.5*d1);
    float f = max(d2/d1,d3/d1,d4/d1);
    printf("%30s ... %s (%3.2fx)\n",name,ok ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m",f);
   	fflush(stdout);
  }
  
  void rtcore_inf(const char* name, RTCSceneFlags sflags, RTCGeometryFlags gflags, int N)
  {
    size_t count = 1000/N;
    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom = addSphere(scene,gflags,zero,2.0f,100);
    rtcCommit (scene);
    size_t numFailures = 0;
    //size_t c0 = __rdtsc();
    double c0 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c1 = __rdtsc();
    double c1 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(inf),dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c2 = __rdtsc();
    double c2 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir+Vec3fa(inf)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c3 = __rdtsc();
    double c3 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(inf),dir+Vec3fa(inf)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c4 = __rdtsc();
    double c4 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir,-0.0f,inf); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    //size_t c5 = __rdtsc();
    double c5 = getSeconds();

    double d1 = c1-c0;
    double d2 = c2-c1;
    double d3 = c3-c2;
    double d4 = c4-c3;
    double d5 = c5-c4;

    rtcDeleteScene (scene);

    bool ok = (d2 < 2.5*d1) && (d3 < 2.5*d1) && (d4 < 2.5*d1) && (d5 < 2.5*d1);
    float f = max(d2/d1,d3/d1,d4/d1,d5/d1);
    printf("%30s ... %s (%3.2fx)\n",name,ok ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m",f);
	  fflush(stdout);
  }

  bool rtcore_overlapping(size_t numTriangles)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, numTriangles, 3);
    AssertNoError();

    Vertex* vertices = (Vertex*) rtcMapBuffer(scene,0,RTC_VERTEX_BUFFER);
    vertices[0].x = 0.0f; vertices[0].y = 0.0f; vertices[0].z = 0.0f;
    vertices[1].x = 1.0f; vertices[1].y = 0.0f; vertices[1].z = 0.0f;
    vertices[2].x = 0.0f; vertices[2].y = 1.0f; vertices[2].z = 0.0f;
    rtcUnmapBuffer(scene,0,RTC_VERTEX_BUFFER);
    AssertNoError();

    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,0,RTC_INDEX_BUFFER);
    for (size_t i=0; i<numTriangles; i++) {
      triangles[i].v0 = 0;
      triangles[i].v1 = 1;
      triangles[i].v2 = 2;
    }
    rtcUnmapBuffer(scene,0,RTC_INDEX_BUFFER);
    AssertNoError();

    rtcCommit (scene);
    AssertNoError();

    return true;
  }

  bool rtcore_backface_culling (RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    /* create triangle that is front facing for a right handed 
     coordinate system if looking along the z direction */
    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned mesh = rtcNewTriangleMesh (scene, gflags, 1, 3);
    Vertex*   vertices  = (Vertex*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    vertices[0].x = 0; vertices[0].y = 0; vertices[0].z = 0;
    vertices[1].x = 0; vertices[1].y = 1; vertices[1].z = 0;
    vertices[2].x = 1; vertices[2].y = 0; vertices[2].z = 0;
    triangles[0].v0 = 0; triangles[0].v1 = 1; triangles[0].v2 = 2;
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    rtcCommit (scene);

    bool passed = true;
    RTCRay ray;
    RTCRay backfacing = makeRay(Vec3fa(0.25f,0.25f,1),Vec3fa(0,0,-1)); 
    RTCRay frontfacing = makeRay(Vec3fa(0.25f,0.25f,-1),Vec3fa(0,0,1)); 

    ray = frontfacing; rtcOccludedN(scene,ray,1);  if (ray.geomID != 0) passed = false;
    ray = frontfacing; rtcIntersectN(scene,ray,1); if (ray.geomID != 0) passed = false;
    ray = backfacing;  rtcOccludedN(scene,ray,1);  if (ray.geomID != -1) passed = false;
    ray = backfacing;  rtcIntersectN(scene,ray,1); if (ray.geomID != -1) passed = false;
#if !defined(__MIC__)
    ray = frontfacing; rtcOccludedN(scene,ray,4);  if (ray.geomID != 0) passed = false;
    ray = frontfacing; rtcIntersectN(scene,ray,4); if (ray.geomID != 0) passed = false;
    ray = backfacing;  rtcOccludedN(scene,ray,4);  if (ray.geomID != -1) passed = false;
    ray = backfacing;  rtcIntersectN(scene,ray,4); if (ray.geomID != -1) passed = false;
#endif
#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    if (has_feature(AVX)) {
      ray = frontfacing; rtcOccludedN(scene,ray,8);  if (ray.geomID != 0) passed = false;
      ray = frontfacing; rtcIntersectN(scene,ray,8); if (ray.geomID != 0) passed = false;
      ray = backfacing;  rtcOccludedN(scene,ray,8);  if (ray.geomID != -1) passed = false;
      ray = backfacing;  rtcIntersectN(scene,ray,8); if (ray.geomID != -1) passed = false;
    }
#endif
#if defined(__MIC__)
    ray = frontfacing; rtcOccludedN(scene,ray,16); if (ray.geomID != 0) passed = false;
    ray = frontfacing; rtcIntersectN(scene,ray,16);if (ray.geomID != 0) passed = false;
    ray = backfacing;  rtcOccludedN(scene,ray,16); if (ray.geomID != -1) passed = false;
    ray = backfacing;  rtcIntersectN(scene,ray,16);if (ray.geomID != -1) passed = false;
#endif
    return passed;
  }

  void rtcore_backface_culling_all ()
  {
    printf("%30s ... ","backface_culling");
    bool passed = true;
    for (int i=0; i<numSceneFlags; i++) 
    {
      RTCSceneFlags flag = getSceneFlag(i);
      bool ok0 = rtcore_backface_culling(flag,RTC_GEOMETRY_STATIC);
      if (ok0) printf("\033[32m+\033[0m"); else printf("\033[31m-\033[0m");
      passed &= ok0;
    }
    printf(" %s\n",passed ? "\033[32m[PASSED]\033[0m" : "\033[31m[FAILED]\033[0m");
	fflush(stdout);
  }

  bool rtcore_new_delete_geometry()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    int geom[1024];
    for (size_t i=0; i<1024; i++) 
      geom[i] = -1;

    for (size_t i=0; i<200; i++) {
      for (size_t j=0; j<20; j++) {
        int index = rand()%1024;
        Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
        if (geom[index] == -1) {
          geom[index] = addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,10);
          AssertNoError();
        }
        else { 
          rtcDeleteGeometry(scene,geom[index]);     
          AssertNoError();
          geom[index] = -1; 
        }
      }
      rtcCommit(scene);
      AssertNoError();
      rtcCommit(scene);
      AssertNoError();
    }
    rtcCommit (scene);
    AssertNoError();
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  void shootRays (RTCScene scene)
  {
    Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
    Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
    RTCRay ray = makeRay(org,dir); 
    rtcOccluded(scene,ray);
    rtcIntersect(scene,ray);

#if !defined(__MIC__)
    RTCRay4 ray4;
    for (size_t j=0; j<4; j++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      setRay(ray4,j,ray);
    }
    __aligned(16) int valid4[4] = { -1,-1,-1,-1 };
    rtcOccluded4(valid4,scene,ray4);
    rtcIntersect4(valid4,scene,ray4);
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    if (has_feature(AVX)) {
      RTCRay8 ray8;
      for (size_t j=0; j<8; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(org,dir); 
        setRay(ray8,j,ray);
      }
      __aligned(32) int valid8[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcOccluded8(valid8,scene,ray8);
      rtcIntersect8(valid8,scene,ray8);
    }
#endif

#if defined(__MIC__)
    RTCRay16 ray16;
    for (size_t j=0; j<16; j++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      setRay(ray16,j,ray);
    }
    __aligned(16) int valid16[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
    rtcOccluded16(valid16,scene,ray16);
    rtcIntersect16(valid16,scene,ray16);
#endif
  }

  bool rtcore_regression_static()
  {
    for (size_t i=0; i<200; i++) 
    {
      if (i%20 == 0) std::cout << "." << std::flush;

      RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
      AssertNoError();

      for (size_t j=0; j<20; j++) 
      {
        Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
        size_t numPhi = rand()%100;
        size_t numTriangles = 2*2*numPhi*(numPhi-1);
        numTriangles = rand()%(numTriangles+1);
        switch (rand()%3) {
        case 0: addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,0.0f); break;
        case 1: addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
        case 2: addUserGeometryEmpty(scene,pos,2.0f); break;
        }
        AssertNoError();
      }

      rtcCommit(scene);
      AssertNoError();

      for (size_t i=0; i<100; i++)
        shootRays(scene);

      rtcDeleteScene (scene);
      AssertNoError();
    }
    return true;
  }

  bool rtcore_regression_dynamic()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    int geom[1024];
    int types[1024];
    size_t numVertices[1024];
    for (size_t i=0; i<1024; i++)  {
      geom[i] = -1;
      types[i] = 0;
      numVertices[i] = 0;
    }

    for (size_t i=0; i<200; i++) 
    {
      if (i%20 == 0) std::cout << "." << std::flush;

      for (size_t j=0; j<20; j++) 
      {
        int index = rand()%1024;
        Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
        if (geom[index] == -1) 
        {
          int type = rand()%3;
          size_t numPhi = rand()%100;
          size_t numTriangles = 2*2*numPhi*(numPhi-1);
          numTriangles = rand()%(numTriangles+1);
          types[index] = type;
          numVertices[index] = 2*numPhi*(numPhi+1);
          switch (type) {
          case 0: geom[index] = addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,0.0f); break;
          case 1: geom[index] = addSphere(scene,RTC_GEOMETRY_DEFORMABLE,pos,2.0f,numPhi,numTriangles,0.0f); break;
          case 2: geom[index] = addSphere(scene,RTC_GEOMETRY_DYNAMIC,pos,2.0f,numPhi,numTriangles,0.0f); break;

          case 3: geom[index] = addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
          case 4: geom[index] = addSphere(scene,RTC_GEOMETRY_DEFORMABLE,pos,2.0f,numPhi,numTriangles,1.0f); break;
          case 5: geom[index] = addSphere(scene,RTC_GEOMETRY_DYNAMIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
            
          case 6: geom[index] = addUserGeometryEmpty(scene,pos,2.0f); break;
          }; 
          AssertNoError();
        }
        else 
        {
          switch (types[index]) {
          case 0:
          case 3:
          case 6: {
            rtcDeleteGeometry(scene,geom[index]);     
            AssertNoError();
            geom[index] = -1; 
            break;
          }
          case 1: 
          case 2:
          case 4: 
          case 5: {
            int op = rand()%2;
            switch (op) {
            case 0: {
              rtcDeleteGeometry(scene,geom[index]);     
              AssertNoError();
              geom[index] = -1; 
              break;
            }
            case 1: {
              Vec3fa* vertices = (Vec3fa*) rtcMapBuffer(scene,geom[index],RTC_VERTEX_BUFFER);
              for (size_t i=0; i<numVertices[index]; i++) vertices[i] += Vec3fa(0.1f);
              rtcUnmapBuffer(scene,geom[index],RTC_VERTEX_BUFFER);

              if (types[index] == 4 || types[index] == 5) {
                Vec3fa* vertices = (Vec3fa*) rtcMapBuffer(scene,geom[index],RTC_VERTEX_BUFFER1);
                for (size_t i=0; i<numVertices[index]; i++) vertices[i] += Vec3fa(0.1f);
                rtcUnmapBuffer(scene,geom[index],RTC_VERTEX_BUFFER1);
              }
              break;
            }
            }
            break;
          }
          }
        }
      }

      rtcCommit(scene);
      AssertNoError();

      for (size_t i=0; i<100; i++)
        shootRays(scene);
    }
    rtcDeleteScene (scene);
    AssertNoError();
    return true;
  }

  /* main function in embree namespace */
  int main(int argc, char** argv) 
  {
    /* parse command line */  
    parseCommandLine(argc,argv);

    /* print Embree version */
    rtcInit("verbose=1");
    rtcExit();

    /* perform tests */
    rtcInit(g_rtcore.c_str());


    POSITIVE("mutex_sys",                 test_mutex_sys());
#if !defined(__MIC__)  // FIXME: hangs on MIC 
    POSITIVE("barrier_sys",               test_barrier_sys());
#endif
#if !defined(__MIC__) && !defined(_WIN32) // FIXME: hangs on MIC and Windows
    POSITIVE("condition_sys",             test_condition_sys());
#endif

    POSITIVE("empty_static",              rtcore_empty(RTC_SCENE_STATIC));
    POSITIVE("empty_dynamic",             rtcore_empty(RTC_SCENE_DYNAMIC));
    POSITIVE("flags_static_static",       rtcore_dynamic_flag(RTC_SCENE_STATIC, RTC_GEOMETRY_STATIC));
    NEGATIVE("flags_static_deformable",   rtcore_dynamic_flag(RTC_SCENE_STATIC, RTC_GEOMETRY_DEFORMABLE));
    NEGATIVE("flags_static_dynamic",      rtcore_dynamic_flag(RTC_SCENE_STATIC, RTC_GEOMETRY_DYNAMIC));
    POSITIVE("flags_dynamic_static",      rtcore_dynamic_flag(RTC_SCENE_DYNAMIC,RTC_GEOMETRY_STATIC));
    POSITIVE("flags_dynamic_deformable",  rtcore_dynamic_flag(RTC_SCENE_DYNAMIC,RTC_GEOMETRY_DEFORMABLE));
    POSITIVE("flags_dynamic_dynamic",     rtcore_dynamic_flag(RTC_SCENE_DYNAMIC,RTC_GEOMETRY_DYNAMIC));
    POSITIVE("static_scene",              rtcore_static_scene());
    //POSITIVE("deformable_geometry",       rtcore_deformable_geometry()); // FIXME
    POSITIVE("unmapped_before_commit",    rtcore_unmapped_before_commit());

#if defined(__BUFFER_STRIDE__)
    POSITIVE("buffer_stride",             rtcore_buffer_stride());
#endif

    POSITIVE("dynamic_enable_disable",    rtcore_dynamic_enable_disable());
    POSITIVE("update_deformable",         rtcore_update(RTC_GEOMETRY_DEFORMABLE));
    POSITIVE("update_dynamic",            rtcore_update(RTC_GEOMETRY_DYNAMIC));
    POSITIVE("overlapping_geometry",      rtcore_overlapping(100000));
    POSITIVE("new_delete_geometry",       rtcore_new_delete_geometry());

#if defined(__USE_RAY_MASK__)
    rtcore_ray_masks_all();
#endif

#if defined(__INTERSECTION_FILTER__)
    rtcore_filter_all();
#endif

#if defined(__BACKFACE_CULLING__)
    rtcore_backface_culling_all();
#endif

    rtcore_packet_write_test_all();

    rtcore_watertight_sphere1(100000);
    rtcore_watertight_plane1(100000);
#if !defined(__MIC__)
    rtcore_watertight_sphere4(100000);
    rtcore_watertight_plane4(100000);
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    if (has_feature(AVX)) {
      rtcore_watertight_sphere8(100000);
      rtcore_watertight_plane8(100000);
    }
#endif

#if defined(__MIC__)
    rtcore_watertight_sphere16(100000);
    rtcore_watertight_plane16(100000);
#endif

#if defined(__FIX_RAYS__)
    rtcore_nan("nan_test_1",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,1);
    rtcore_inf("inf_test_1",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,1);

#if !defined(__MIC__)
    rtcore_nan("nan_test_4",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,4);
    rtcore_inf("inf_test_4",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,4);
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    if (has_feature(AVX)) {
      rtcore_nan("nan_test_8",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,8);
      rtcore_inf("inf_test_8",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,8);
    }
#endif

#if defined(__MIC__)
    rtcore_nan("nan_test_16",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,16);
    rtcore_inf("inf_test_16",RTC_SCENE_STATIC,RTC_GEOMETRY_STATIC,16);
#endif
#endif

    POSITIVE("regression_static",         rtcore_regression_static());
    POSITIVE("regression_dynamic",        rtcore_regression_dynamic());

    rtcExit();

    return 0;
  }
}

int main(int argc, char** argv)
{
  try {
    return embree::main(argc, argv);
  }
  catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Error: unknown exception caught." << std::endl;
    return 1;
  }
}

