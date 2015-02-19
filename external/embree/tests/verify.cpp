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

//#define DEFAULT_STACK_SIZE 2*1024*1024
//#define DEFAULT_STACK_SIZE 512*1024
#define DEFAULT_STACK_SIZE 0

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
  //static size_t testN = 10000000;
  static size_t regressionN = 200;

  /* vertex and triangle layout */
  struct Vertex  { float x,y,z,a; };
#if defined(__MIC__)
  typedef Vec3fa Vertex3f;
  typedef Vec3fa Vertex3fa;
#else
  typedef Vec3f  Vertex3f;
  typedef Vec3fa Vertex3fa;
#endif
  struct Triangle { int v0, v1, v2; };

  std::vector<thread_t> g_threads;
  size_t numFailedTests = 0;

  atomic_t errorCounter = 0;

#if defined(__WIN32__)
#  define GREEN(x) x
#  define RED(x) x
#else
#  define GREEN(x) "\033[32m" x "\033[0m"
#  define RED(x) "\033[31m" x "\033[0m"
#endif

#define CountErrors() \
  if (rtcGetError() != RTC_NO_ERROR) atomic_add(&errorCounter,1);
#define AssertNoError() \
  if (rtcGetError() != RTC_NO_ERROR) return false;
#define AssertAnyError() \
  if (rtcGetError() == RTC_NO_ERROR) return false;
#define AssertError(code) \
  if (rtcGetError() != code) return false;
#define POSITIVE(name,test) {                                               \
    printf("%30s ...",name);                                            \
    bool ok = test;                                                     \
    printf(" %s\n",ok ? GREEN("[PASSED]") : RED("[FAILED]"));          \
    fflush(stdout);                                                     \
    numFailedTests += !ok;                                              \
  }
#define NEGATIVE(name,test) {                                       \
    printf("%30s ... ",name);                                           \
    bool notok = test;                                                  \
    printf(" %s\n",notok ? RED("[FAILED]") : GREEN("[PASSED]")); \
    fflush(stdout);                                                     \
    numFailedTests += notok;                                              \
  }

  const size_t numSceneFlags = 64;

  std::vector<void*> buffers;
  MutexSys g_mutex2;
  void* allocBuffer(size_t size) { 
    g_mutex2.lock();
    void* ptr = alignedMalloc(size);
    buffers.push_back(ptr); 
    g_mutex2.unlock();
    return ptr; 
  }
  void clearBuffers() {
    for (size_t i=0; i<buffers.size(); i++) {
      alignedFree(buffers[i]);
    }
    buffers.clear();
  }

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

  const size_t numSceneGeomFlags = 32;

  void getSceneGeomFlag(size_t i, RTCSceneFlags& sflags, RTCGeometryFlags& gflags) 
  {
    int sflag = 0, gflag = 0;
    if (i & 4) {
      sflag |= RTC_SCENE_DYNAMIC;
      gflag = min(i&3,size_t(2));
    }
    if (i & 8) sflag |= RTC_SCENE_HIGH_QUALITY;
    if (i & 16) sflag |= RTC_SCENE_ROBUST;
    sflags = (RTCSceneFlags) sflag;
    gflags = (RTCGeometryFlags) gflag;
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

      /* get number of regression test iterations to perform */
      else if (tag == "-regressions") {
        if (i+1 >= argc) THROW_RUNTIME_ERROR("command line parsing error");
        regressionN = atoi(argv[++i]);
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
    Vertex3fa*   vertices  = (Vertex3fa*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
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
    size_t numVertices = numTheta*(numPhi+1);
    
    unsigned mesh = rtcNewTriangleMesh (scene, flag, numTriangles, numVertices,numTimeSteps);
    
    /* map triangle and vertex buffer */
    Vertex3f* vertices0 = NULL;
    Vertex3f* vertices1 = NULL;
    if (numTimeSteps >= 1) rtcSetBuffer(scene,mesh,RTC_VERTEX_BUFFER0,vertices0 = (Vertex3f*) allocBuffer(numVertices*sizeof(Vertex3f)), 0, sizeof(Vertex3f)); 
    if (numTimeSteps >= 2) rtcSetBuffer(scene,mesh,RTC_VERTEX_BUFFER1,vertices1 = (Vertex3f*) allocBuffer(numVertices*sizeof(Vertex3f)), 0, sizeof(Vertex3f)); 
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
        Vertex3f* v = &vertices0[phi*numTheta+theta];
        const float cosThetaf = cos(thetaf);
        v->x = pos.x + r*sin(phif)*sin(thetaf);
        v->y = pos.y + r*cos(phif);
        v->z = pos.z + r*sin(phif)*cosThetaf;

        if (vertices1) {
          Vertex3f* v1 = &vertices1[phi*numTheta+theta];
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

    //if (numTimeSteps >= 1) rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
    //if (numTimeSteps >= 2) rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    return mesh;
  }

  /* adds a subdiv sphere to the scene */
  unsigned int addSubdivSphere (RTCScene scene, RTCGeometryFlags flags, const Vec3fa& pos, const float r, size_t numPhi, float level, size_t maxFaces = -1, float motion = 0.0f)
  {
    size_t numTheta = 2*numPhi;
    vector_t<Vec3fa> vertices(numTheta*(numPhi+1));
    std::vector<int> indices;
    std::vector<int> faces;
    std::vector<int> offsets;
    
    /* create sphere geometry */
    const float rcpNumTheta = rcp((float)numTheta);
    const float rcpNumPhi   = rcp((float)numPhi);
    for (int phi=0; phi<=numPhi; phi++)
    {
      for (int theta=0; theta<numTheta; theta++)
      {
	const float phif   = phi*float(pi)*rcpNumPhi;
	const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
	Vec3fa& v = vertices[phi*numTheta+theta];
	Vec3fa P(pos.x + r*sin(phif)*sin(thetaf),
		 pos.y + r*cos(phif),
		 pos.z + r*sin(phif)*cos(thetaf));
	v.x = P.x;
	v.y = P.y;
	v.z = P.z;
      }
      if (phi == 0) continue;
      
      if (phi == 1)
      {
	for (int theta=1; theta<=numTheta; theta++) 
	{
	  int p00 = numTheta-1;
	  int p10 = phi*numTheta+theta-1;
	  int p11 = phi*numTheta+theta%numTheta;
	  offsets.push_back(indices.size());
	  indices.push_back(p10); 
	  indices.push_back(p00);
	  indices.push_back(p11);
	  faces.push_back(3);
	}
      }
      else if (phi == numPhi)
      {
	for (int theta=1; theta<=numTheta; theta++) 
	{
	  int p00 = (phi-1)*numTheta+theta-1;
	  int p01 = (phi-1)*numTheta+theta%numTheta;
	  int p10 = numPhi*numTheta;
	  offsets.push_back(indices.size());
	  indices.push_back(p10);
	  indices.push_back(p00);
	  indices.push_back(p01);
	  faces.push_back(3);
	}
      }
      else
      {
	for (int theta=1; theta<=numTheta; theta++) 
	{
	  int p00 = (phi-1)*numTheta+theta-1;
	  int p01 = (phi-1)*numTheta+theta%numTheta;
	  int p10 = phi*numTheta+theta-1;
	  int p11 = phi*numTheta+theta%numTheta;
	  offsets.push_back(indices.size());
	  indices.push_back(p10);
	  indices.push_back(p00);
	  indices.push_back(p01);
	  indices.push_back(p11);
	  faces.push_back(4);
	}
      }
    }
    
    /* create subdiv geometry */
    size_t numFaces = min(faces.size(),maxFaces);
    size_t numEdges = indices.size();
    size_t numVertices = vertices.size();
    size_t numEdgeCreases = 10;
    size_t numVertexCreases = 10;
    size_t numHoles = 0; // do not test holes as this causes some tests that assume a closed sphere to fail
    unsigned int mesh = rtcNewSubdivisionMesh(scene, flags, numFaces, numEdges, numVertices, numEdgeCreases, numVertexCreases, numHoles);
    Vec3fa* vertexBuffer = (Vec3fa*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    int*    indexBuffer  = (int     *) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    int*    facesBuffer = (int     *) rtcMapBuffer(scene,mesh,RTC_FACE_BUFFER);
    float*  levelBuffer  = (float   *) rtcMapBuffer(scene,mesh,RTC_LEVEL_BUFFER);

    memcpy(vertexBuffer,vertices.data(),numVertices*sizeof(Vec3fa));
    memcpy(indexBuffer ,indices.data() ,numEdges*sizeof(int));
    memcpy(facesBuffer,faces.data() ,numFaces*sizeof(int));
    for (size_t i=0; i<indices.size(); i++) levelBuffer[i] = level;
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    rtcUnmapBuffer(scene,mesh,RTC_FACE_BUFFER);
    rtcUnmapBuffer(scene,mesh,RTC_LEVEL_BUFFER);
    
    int* edgeCreaseIndices  = (int*) rtcMapBuffer(scene,mesh,RTC_EDGE_CREASE_INDEX_BUFFER);
    float* edgeCreaseWeights = (float*) rtcMapBuffer(scene,mesh,RTC_EDGE_CREASE_WEIGHT_BUFFER);
    for (size_t i=0; i<numEdgeCreases; i++) 
    {
      if (faces.size()) {
	int f = rand() % faces.size();
	int n = faces[f];
	int e = rand() % n;
	edgeCreaseIndices[2*i+0] = indices[offsets[f]+(e+0)%n];
	edgeCreaseIndices[2*i+1] = indices[offsets[f]+(e+1)%n];
      } else {
	edgeCreaseIndices[2*i+0] = 0;
	edgeCreaseIndices[2*i+1] = 0;
      }
      edgeCreaseWeights[i] = 10.0f*drand48();
    }
    rtcUnmapBuffer(scene,mesh,RTC_EDGE_CREASE_INDEX_BUFFER); 
    rtcUnmapBuffer(scene,mesh,RTC_EDGE_CREASE_WEIGHT_BUFFER); 
    
    int* vertexCreaseIndices  = (int*) rtcMapBuffer(scene,mesh,RTC_VERTEX_CREASE_INDEX_BUFFER);
    float* vertexCreaseWeights = (float*) rtcMapBuffer(scene,mesh,RTC_VERTEX_CREASE_WEIGHT_BUFFER);
    for (size_t i=0; i<numVertexCreases; i++) 
    {
      int v = numTheta-1 + rand() % (vertices.size()+2-2*numTheta);
      vertexCreaseIndices[i] = v;
      vertexCreaseWeights[i] = 10.0f*drand48();
    }
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_CREASE_INDEX_BUFFER); 
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_CREASE_WEIGHT_BUFFER); 
    
    int* holeBuffer  = (int*) rtcMapBuffer(scene,mesh,RTC_HOLE_BUFFER);
    for (size_t i=0; i<numHoles; i++) {
      holeBuffer[i] = rand() % faces.size();
    }
    rtcUnmapBuffer(scene,mesh,RTC_HOLE_BUFFER); 
    
    return mesh;
  }

  unsigned int addCube (RTCScene scene_i, RTCGeometryFlags flag, const Vec3fa& pos, const float r)
  {
    /* create a triangulated cube with 12 triangles and 8 vertices */
    unsigned int mesh = rtcNewTriangleMesh (scene_i, flag, 12, 8);
    
    /* set vertices */
    Vec3fa* vertices = (Vec3fa*) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
    vertices[0] = pos + r*Vec3fa(-1,-1,-1); 
    vertices[1] = pos + r*Vec3fa(-1,-1,+1); 
    vertices[2] = pos + r*Vec3fa(-1,+1,-1); 
    vertices[3] = pos + r*Vec3fa(-1,+1,+1); 
    vertices[4] = pos + r*Vec3fa(+1,-1,-1); 
    vertices[5] = pos + r*Vec3fa(+1,-1,+1); 
    vertices[6] = pos + r*Vec3fa(+1,+1,-1); 
    vertices[7] = pos + r*Vec3fa(+1,+1,+1); 
    rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

    /* set triangles and colors */
    int tri = 0;
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
    
    // left side
    triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; tri++;
    triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; tri++;

    // right side
    triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; tri++;
    triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;
    
    // bottom side
    triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; tri++;
    triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; tri++;
    
    // top side
    triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; tri++;
    triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;
    
    // front side
    triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; tri++;
    triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; tri++;
    
    // back side
    triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; tri++;
    triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; tri++;
    
    rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
    
    return mesh;
  }

  unsigned addHair (RTCScene scene, RTCGeometryFlags flag, const Vec3fa& pos, const float scale, const float r, size_t numHairs = 1, float motion = 0.0f)
  {
    size_t numTimeSteps = motion == 0.0f ? 1 : 2;
    unsigned geomID = rtcNewHairGeometry (scene, flag, numHairs, numHairs*4, numTimeSteps);
    
    /* map triangle and vertex buffer */
    Vec3fa* vertices0 = NULL;
    Vec3fa* vertices1 = NULL;
    if (numTimeSteps >= 1) vertices0 = (Vec3fa*) rtcMapBuffer(scene,geomID,RTC_VERTEX_BUFFER0); 
    if (numTimeSteps >= 2) vertices1 = (Vec3fa*) rtcMapBuffer(scene,geomID,RTC_VERTEX_BUFFER1); 
    int* indices = (int*) rtcMapBuffer(scene,geomID,RTC_INDEX_BUFFER);

    for (size_t i=0; i<numHairs; i++) 
    {
      indices[i] = 4*i;
      const Vec3fa p0 = pos + scale*Vec3fa(i%7,i%13,i%31);
      const Vec3fa p1 = p0 + scale*Vec3fa(1,0,0);
      const Vec3fa p2 = p0 + scale*Vec3fa(0,1,1);
      const Vec3fa p3 = p0 + scale*Vec3fa(0,1,0);
      
      if (vertices0) {
        vertices0[4*i+0] = Vec3fa(p0,r);
        vertices0[4*i+1] = Vec3fa(p1,r);
        vertices0[4*i+2] = Vec3fa(p2,r);
        vertices0[4*i+3] = Vec3fa(p3,r);
      }
      if (vertices1) {
        vertices1[4*i+0] = Vec3fa(p0+Vec3fa(motion),r);
        vertices1[4*i+1] = Vec3fa(p1+Vec3fa(motion),r);
        vertices1[4*i+2] = Vec3fa(p2+Vec3fa(motion),r);
        vertices1[4*i+3] = Vec3fa(p3+Vec3fa(motion),r);
      }
    }

    if (numTimeSteps >= 1) rtcUnmapBuffer(scene,geomID,RTC_VERTEX_BUFFER0); 
    if (numTimeSteps >= 2) rtcUnmapBuffer(scene,geomID,RTC_VERTEX_BUFFER1); 
    rtcUnmapBuffer(scene,geomID,RTC_INDEX_BUFFER);
    return geomID;
  }

  unsigned addGarbageTriangles (RTCScene scene, RTCGeometryFlags flag, size_t numTriangles, bool motion)
  {
    /* create a triangulated sphere */
    size_t numTimeSteps = motion ? 1 : 2;
    unsigned mesh = rtcNewTriangleMesh (scene, flag, numTriangles, 3*numTriangles,numTimeSteps);
    
    /* map triangle and vertex buffer */
    if (numTimeSteps >= 1) {
      int* v = (int*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
      for (size_t i=0; i<4*3*numTriangles; i++) v[i] = rand();
      rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
    }
    if (numTimeSteps >= 2) {
      int* v = (int*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
      for (size_t i=0; i<4*3*numTriangles; i++) v[i] = rand();
      rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
    }
    
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    for (size_t i=0; i<numTriangles; i++) {
      triangles[i].v0 = (rand() % 32 == 0) ? rand() : 3*i+0;
      triangles[i].v1 = (rand() % 32 == 0) ? rand() : 3*i+1;
      triangles[i].v2 = (rand() % 32 == 0) ? rand() : 3*i+2;
    }
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);

    return mesh;
  }

  unsigned addGarbageHair (RTCScene scene, RTCGeometryFlags flag, size_t numCurves, bool motion)
  {
    /* create a triangulated sphere */
    size_t numTimeSteps = motion ? 1 : 2;
    unsigned mesh = rtcNewHairGeometry (scene, flag, numCurves, 4*numCurves,numTimeSteps);
    
    /* map triangle and vertex buffer */
    if (numTimeSteps >= 1) {
      int* v = (int*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
      for (size_t i=0; i<4*4*numCurves; i++) v[i] = rand();
      rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER0); 
    }
    if (numTimeSteps >= 2) {
      int* v = (int*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
      for (size_t i=0; i<4*4*numCurves; i++) v[i] = rand();
      rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER1); 
    }
    
    int* curves = (int*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
    for (size_t i=0; i<numCurves; i++) 
      curves[i] = (rand() % 32 == 0) ? rand() : 4*i;
    rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);

    return mesh;
  }

  struct Sphere
  {
    ALIGNED_CLASS;
  public:
    Sphere () : pos(zero), r(zero) {}
    Sphere (const Vec3fa& pos, float r) : pos(pos), r(r) {}
    __forceinline BBox3fa bounds() const { return BBox3fa(pos-Vec3fa(r),pos+Vec3fa(r)); }
  public:
    Vec3fa pos;
    float r;
  };

  void BoundsFunc(Sphere* sphere, size_t index, BBox3fa* bounds_o)
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

  unsigned addUserGeometryEmpty (RTCScene scene, Sphere* sphere)
  {
    BBox3fa bounds = sphere->bounds(); 
    unsigned geom = rtcNewUserGeometry (scene,1);
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
    for (size_t i=0; i<1000; i++) 
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
      g_threads.push_back(createThread(test_barrier_sys_thread,NULL,DEFAULT_STACK_SIZE,i));
    setAffinity(0);
    
    bool ok = true;
    for (size_t i=0; i<1000; i++) 
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
      g_threads.push_back(createThread(test_condition_sys_thread,NULL,DEFAULT_STACK_SIZE,i));
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
    numFailedTests += !ok;
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
      g_threads.push_back(createThread(test_mutex_sys_thread,NULL,DEFAULT_STACK_SIZE,i));

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
    bool ok = g_counter == 10000*numThreads;
    numFailedTests += !ok;
    return ok;
  }

  bool rtcore_empty(RTCSceneFlags flags)
  {
    RTCScene scene = rtcNewScene(flags,aflags);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
    rtcDeleteScene (scene);
    clearBuffers();
    AssertNoError();
    return true;
  }

  bool rtcore_dynamic_flag(RTCSceneFlags sceneFlag, RTCGeometryFlags geomFlag)
  {
    RTCScene scene = rtcNewScene(sceneFlag,aflags);
    AssertNoError();
    rtcNewTriangleMesh (scene, geomFlag, 0, 0);
    AssertNoError();
    rtcNewHairGeometry (scene, geomFlag, 0, 0);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
    rtcDeleteScene (scene);
    clearBuffers();
    AssertNoError();
    return true;
  }

  bool rtcore_static_scene()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    unsigned geom0 = addSphere(scene,RTC_GEOMETRY_STATIC,zero,1.0f,50);
    AssertNoError();
    unsigned geom1 = addSubdivSphere(scene,RTC_GEOMETRY_STATIC,zero,1.0f,10,16);
    AssertNoError();
    unsigned geom2 = addHair(scene,RTC_GEOMETRY_STATIC,Vec3fa(0,0,0),1.0f,0.5f,100);
    AssertNoError();
    rtcCommit (scene);
    AssertNoError();
    rtcCommit (scene); // cannot commit static scene twice
    AssertAnyError();
    rtcDisable(scene,geom0); // static scene cannot get modified anymore after commit
    AssertAnyError();
    rtcDeleteScene (scene);
    clearBuffers();
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
    rtcMapBuffer(scene,geom,RTC_INDEX_BUFFER);
    AssertError(RTC_INVALID_OPERATION); // cannot modify index buffer of deformable geometry anymore after commit
    rtcMapBuffer(scene,geom,RTC_VERTEX_BUFFER);
    AssertNoError();
    rtcUnmapBuffer(scene,geom,RTC_INDEX_BUFFER);
    AssertError(RTC_INVALID_OPERATION); // cannot modify index buffer of deformable geometry anymore after commit
    rtcUnmapBuffer(scene,geom,RTC_VERTEX_BUFFER);
    AssertNoError();
    rtcDeleteScene (scene);
    clearBuffers();
    AssertNoError();
    return true;
  }

  bool rtcore_unmapped_before_commit()
  {
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
    clearBuffers();
    AssertNoError();
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
    clearBuffers();
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
    //unsigned geom1 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(-1,0,+1),1.0f,50);
    unsigned geom1 = addHair  (scene,RTC_GEOMETRY_STATIC,Vec3fa(-1,0,+1),1.0f,1.0f,1);
    unsigned geom2 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(+1,0,-1),1.0f,50);
    //unsigned geom3 = addSphere(scene,RTC_GEOMETRY_STATIC,Vec3fa(+1,0,+1),1.0f,50);
    unsigned geom3 = addHair  (scene,RTC_GEOMETRY_STATIC,Vec3fa(+1,0,+1),1.0f,1.0f,1);
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
    clearBuffers();
    AssertNoError();
    return true;
  }

  void move_mesh_vec3f(RTCScene scene, unsigned mesh, size_t numVertices, Vec3fa& pos) 
  {
    Vertex3f* vertices = (Vertex3f*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    for (size_t i=0; i<numVertices; i++) vertices[i] += Vertex3f(pos);
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER);
    rtcUpdate(scene,mesh);
  }

  void move_mesh_vec3fa(RTCScene scene, unsigned mesh, size_t numVertices, Vec3fa& pos) 
  {
    Vertex3fa* vertices = (Vertex3fa*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
    for (size_t i=0; i<numVertices; i++) vertices[i] += Vertex3fa(pos);
    rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER);
    rtcUpdate(scene,mesh);
  }
  
  bool rtcore_update(RTCGeometryFlags flags)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    size_t numPhi = 50;
    size_t numVertices = 2*numPhi*(numPhi+1);
    Vec3fa pos0 = Vec3fa(-10,0,-10);
    Vec3fa pos1 = Vec3fa(-10,0,+10);
    Vec3fa pos2 = Vec3fa(+10,0,-10);
    Vec3fa pos3 = Vec3fa(+10,0,+10);
    unsigned geom0 = addSphere(scene,flags,pos0,1.0f,numPhi);
    unsigned geom1 = addHair  (scene,flags,pos1,1.0f,1.0f,1);
    unsigned geom2 = addSphere(scene,flags,pos2,1.0f,numPhi);
    unsigned geom3 = addHair  (scene,flags,pos3,1.0f,1.0f,1);
    AssertNoError();
    
    for (size_t i=0; i<16; i++) 
    {
      bool move0 = i & 1, move1 = i & 2, move2 = i & 4, move3 = i & 8;
      Vec3fa ds(2,0.1f,2);
      if (move0) { move_mesh_vec3f (scene,geom0,numVertices,ds); pos0 += ds; }
      if (move1) { move_mesh_vec3fa(scene,geom1,4,ds); pos1 += ds; }
      if (move2) { move_mesh_vec3f (scene,geom2,numVertices,ds); pos2 += ds; }
      if (move3) { move_mesh_vec3fa(scene,geom3,4,ds); pos3 += ds; }
      rtcCommit (scene);
      AssertNoError();
      {
        RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); 
        RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); 
        RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); 
        RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); 
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
    clearBuffers();
    AssertNoError();
    return true;
  }

  bool rtcore_ray_masks_intersect(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;
    Vec3fa pos0 = Vec3fa(-10,0,-10);
    Vec3fa pos1 = Vec3fa(-10,0,+10);
    Vec3fa pos2 = Vec3fa(+10,0,-10);
    Vec3fa pos3 = Vec3fa(+10,0,+10);

    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom0 = addSphere(scene,gflags,pos0,1.0f,50);
    //unsigned geom1 = addSphere(scene,gflags,pos1,1.0f,50);
    unsigned geom1 = addHair  (scene,gflags,pos1,1.0f,1.0f,1);
    unsigned geom2 = addSphere(scene,gflags,pos2,1.0f,50);
    //unsigned geom3 = addSphere(scene,gflags,pos3,1.0f,50);
    unsigned geom3 = addHair  (scene,gflags,pos3,1.0f,1.0f,1);
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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;
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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
    clearBuffers();
    return passed;
  }

  bool rtcore_ray_masks_occluded(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    bool passed = true;
    Vec3fa pos0 = Vec3fa(-10,0,-10);
    Vec3fa pos1 = Vec3fa(-10,0,+10);
    Vec3fa pos2 = Vec3fa(+10,0,-10);
    Vec3fa pos3 = Vec3fa(+10,0,+10);

    RTCScene scene = rtcNewScene(sflags,aflags);
    unsigned geom0 = addSphere(scene,gflags,pos0,1.0f,50);
    unsigned geom1 = addSphere(scene,gflags,pos1,1.0f,50);
    unsigned geom2 = addSphere(scene,gflags,pos2,1.0f,50);
    unsigned geom3 = addSphere(scene,gflags,pos3,1.0f,50);
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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;
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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
	RTCRay ray0 = makeRay(pos0+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray0.mask = mask0;
	RTCRay ray1 = makeRay(pos1+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray1.mask = mask1;
	RTCRay ray2 = makeRay(pos2+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray2.mask = mask2;
	RTCRay ray3 = makeRay(pos3+Vec3fa(0,10,0),Vec3fa(0,-1,0)); ray3.mask = mask3;

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
    clearBuffers();
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
      if (ok0) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok0;
      bool ok1 = rtcore_ray_masks_occluded(flag,RTC_GEOMETRY_STATIC);
      if (ok1) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok1;
    }
    printf(" %s\n",passed ? GREEN("[PASSED]") : RED("[FAILED]"));
    fflush(stdout);
    numFailedTests += !passed;
  }

  bool rtcore_build(RTCSceneFlags sflags, RTCGeometryFlags gflags)
  {
    RTCScene scene = rtcNewScene(sflags,aflags);
    addSphere(scene,gflags,zero,1E-24f,50);
    addHair(scene,gflags,zero,1E-24f,1E-26f,100,1E-26f);
    addSphere(scene,gflags,zero,1E-24f,50);
    addHair(scene,gflags,zero,1E-24f,1E-26f,100,1E-26f);
    rtcCommit (scene);
    rtcDeleteScene (scene);
    clearBuffers();
    return true;
  }
  
  void rtcore_build()
  {
    printf("%30s ... ","build");
    bool passed = true;
    for (int i=0; i<numSceneGeomFlags; i++) 
    {
      RTCSceneFlags sflags; RTCGeometryFlags gflags;
      getSceneGeomFlag(i,sflags,gflags);
      bool ok = rtcore_build(sflags,gflags);
      if (ok) printf(GREEN("+")); else printf(RED("-"));
    }
    printf(" %s\n",true ? GREEN("[PASSED]") : RED("[FAILED]"));
    fflush(stdout);
    numFailedTests += !passed;
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

#if defined(__MIC__)
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
    clearBuffers();
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
    clearBuffers();
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
      if (ok0) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok0;
      bool ok1 = rtcore_filter_occluded(flag,RTC_GEOMETRY_STATIC);
      if (ok1) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok1;

    }
    printf(" %s\n",passed ? GREEN("[PASSED]") : RED("[FAILED]"));
    fflush(stdout);
    numFailedTests += !passed;
  }

  bool rtcore_packet_write_test(RTCSceneFlags sflags, RTCGeometryFlags gflags, int type)
  {
    bool passed = true;
    RTCScene scene = rtcNewScene(sflags,aflags);

    switch (type) {
    case 0: addSphere(scene,gflags,Vec3fa(-1,0,-1),1.0f,50,-1,0.0f); break;
    case 1: addSphere(scene,gflags,Vec3fa(-1,0,-1),1.0f,50,-1,0.1f); break;
    case 2: addHair  (scene,gflags,Vec3fa(-1,0,-1),1.0f,1.0f,1,0.0f); break;
    case 3: addHair  (scene,gflags,Vec3fa(-1,0,-1),1.0f,1.0f,1,0.1f); break; 
    }
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
      bool ok = true;
      ok &= rtcore_packet_write_test(flag,RTC_GEOMETRY_STATIC,0);
      ok &= rtcore_packet_write_test(flag,RTC_GEOMETRY_STATIC,1);
      ok &= rtcore_packet_write_test(flag,RTC_GEOMETRY_STATIC,2);
      ok &= rtcore_packet_write_test(flag,RTC_GEOMETRY_STATIC,3);
      if (ok) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok;
    }
    printf(" %s\n",passed ? GREEN("[PASSED]") : RED("[FAILED]"));
    fflush(stdout);
    numFailedTests += !passed;
  }

  void rtcore_watertight_closed1(const std::string& type, const Vec3fa& pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    if      (type == "sphere") addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,500);
    else if (type == "cube"  ) addCube  (scene,RTC_GEOMETRY_STATIC,pos,2.0f);
    rtcCommit (scene);

    size_t numFailures = 0;
    for (size_t i=0; i<testN; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(pos+org,dir); 
      rtcIntersect(scene,ray);
      numFailures += ray.primID == -1;
    }
    rtcDeleteScene (scene);
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n", ("watertight_"+type+"1").c_str(), failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }
  
  void rtcore_watertight_closed4(const std::string& type, const Vec3fa& pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    if      (type == "sphere") addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,500);
    else if (type == "cube"  ) addCube  (scene,RTC_GEOMETRY_STATIC,pos,2.0f);
    rtcCommit (scene);

    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=4) {
      RTCRay4 ray4;
      for (size_t j=0; j<4; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(pos+org,dir); 
        setRay(ray4,j,ray);
      }
      __aligned(16) int valid[4] = { -1,-1,-1,-1 };
      rtcIntersect4(valid,scene,ray4);
      for (size_t j=0; j<4; j++)
        numFailures += ray4.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n", ("watertight_"+type+"4").c_str(), failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_watertight_closed8(const std::string& type, const Vec3fa& pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    if      (type == "sphere") addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,500);
    else if (type == "cube"  ) addCube  (scene,RTC_GEOMETRY_STATIC,pos,2.0f);
    rtcCommit (scene);

    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=8) {
      RTCRay8 ray8;
      for (size_t j=0; j<8; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(pos+org,dir); 
        setRay(ray8,j,ray);
      }
      __aligned(32) int valid[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect8(valid,scene,ray8);
      for (size_t j=0; j<8; j++)
        numFailures += ray8.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n", ("watertight_"+type+"8").c_str(), failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_watertight_closed16(const std::string& type, const Vec3fa& pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    if      (type == "sphere") addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,500);
    else if (type == "cube"  ) addCube  (scene,RTC_GEOMETRY_STATIC,pos,2.0f);
    rtcCommit (scene);

    size_t numFailures = 0;
    for (size_t i=0; i<testN; i+=16) {
      RTCRay16 ray16;
      for (size_t j=0; j<16; j++) {
        Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
        RTCRay ray = makeRay(pos+org,dir); 
        setRay(ray16,j,ray);
      }
      __aligned(64) int valid[16] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
      rtcIntersect16(valid,scene,ray16);
      for (size_t j=0; j<16; j++)
        numFailures += ray16.primID[j] == -1;
    }
    rtcDeleteScene (scene);
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n", ("watertight_"+type+"16").c_str(), failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }
  
  void rtcore_watertight_plane1(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,500,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
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
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n","watertight_plane1", failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_watertight_plane4(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,500,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
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
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n","watertight_plane4", failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_watertight_plane8(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,500,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
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
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n","watertight_plane8",failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_watertight_plane16(float pos)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_ROBUST,aflags);
    unsigned geom = addPlane(scene,RTC_GEOMETRY_STATIC,500,Vec3fa(pos,-6.0f,-6.0f),Vec3fa(0.0f,12.0f,0.0f),Vec3fa(0.0f,0.0f,12.0f));
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
    clearBuffers();
    double failRate = double(numFailures) / double(testN);
    bool failed = failRate > 0.00002;

    printf("%30s ... %s (%f%%)\n","watertight_plane16", failed ? RED("[FAILED]") : GREEN("[PASSED]"), 100.0f*failRate);
    fflush(stdout);
    numFailedTests += failed;
  }

  void rtcore_nan(const char* name, RTCSceneFlags sflags, RTCGeometryFlags gflags, int N)
  {
    size_t count = 1000/N;
    RTCScene scene = rtcNewScene(sflags,aflags);
    addSphere(scene,gflags,zero,2.0f,100);
    addHair  (scene,gflags,zero,1.0f,1.0f,100);
    rtcCommit (scene);
    size_t numFailures = 0;
    double c0 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c1 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(nan),dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c2 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(nan),dir+Vec3fa(nan)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c3 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir,nan,nan); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c4 = getSeconds();

    double d1 = c1-c0;
    double d2 = c2-c1;
    double d3 = c3-c2;
    double d4 = c4-c3;
    rtcDeleteScene (scene);
    clearBuffers();

    bool ok = (d2 < 2.5*d1) && (d3 < 2.5*d1) && (d4 < 2.5*d1);
    float f = max(d2/d1,d3/d1,d4/d1);
    printf("%30s ... %s (%3.2fx)\n",name,ok ? GREEN("[PASSED]") : RED("[FAILED]"),f);
    fflush(stdout);
    numFailedTests += !ok;
  }
  
  void rtcore_inf(const char* name, RTCSceneFlags sflags, RTCGeometryFlags gflags, int N)
  {
    size_t count = 1000/N;
    RTCScene scene = rtcNewScene(sflags,aflags);
    addSphere(scene,gflags,zero,2.0f,100);
    addHair  (scene,gflags,zero,1.0f,1.0f,100);
    rtcCommit (scene);
    size_t numFailures = 0;
    double c0 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c1 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(inf),dir); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c2 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir+Vec3fa(inf)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c3 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org+Vec3fa(inf),dir+Vec3fa(inf)); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c4 = getSeconds();
    for (size_t i=0; i<count; i++) {
      Vec3fa org(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      Vec3fa dir(2.0f*drand48()-1.0f,2.0f*drand48()-1.0f,2.0f*drand48()-1.0f);
      RTCRay ray = makeRay(org,dir,-0.0f,inf); 
      rtcOccludedN(scene,ray,N);
      rtcIntersectN(scene,ray,N);
    }
    double c5 = getSeconds();

    double d1 = c1-c0;
    double d2 = c2-c1;
    double d3 = c3-c2;
    double d4 = c4-c3;
    double d5 = c5-c4;

    rtcDeleteScene (scene);
    clearBuffers();

    bool ok = (d2 < 2.5*d1) && (d3 < 2.5*d1) && (d4 < 2.5*d1) && (d5 < 2.5*d1);
    float f = max(d2/d1,d3/d1,d4/d1,d5/d1);
    printf("%30s ... %s (%3.2fx)\n",name,ok ? GREEN("[PASSED]") : RED("[FAILED]"),f);
    fflush(stdout);
    numFailedTests += !ok;
  }

  bool rtcore_overlapping_triangles(size_t N)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, N, 3);
    AssertNoError();

    Vertex3fa* vertices = (Vertex3fa*) rtcMapBuffer(scene,0,RTC_VERTEX_BUFFER);
    vertices[0].x = 0.0f; vertices[0].y = 0.0f; vertices[0].z = 0.0f;
    vertices[1].x = 1.0f; vertices[1].y = 0.0f; vertices[1].z = 0.0f;
    vertices[2].x = 0.0f; vertices[2].y = 1.0f; vertices[2].z = 0.0f;
    rtcUnmapBuffer(scene,0,RTC_VERTEX_BUFFER);
    AssertNoError();

    Triangle* triangles = (Triangle*) rtcMapBuffer(scene,0,RTC_INDEX_BUFFER);
    for (size_t i=0; i<N; i++) {
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

  bool rtcore_overlapping_hair(size_t N)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);
    AssertNoError();
    rtcNewHairGeometry (scene, RTC_GEOMETRY_STATIC, N, 4);
    AssertNoError();

    Vec3fa* vertices = (Vec3fa*) rtcMapBuffer(scene,0,RTC_VERTEX_BUFFER);
    vertices[0].x = 0.0f; vertices[0].y = 0.0f; vertices[0].z = 0.0f; vertices[0].w = 0.1f;
    vertices[1].x = 0.0f; vertices[1].y = 0.0f; vertices[1].z = 1.0f; vertices[1].w = 0.1f;
    vertices[2].x = 0.0f; vertices[2].y = 1.0f; vertices[2].z = 1.0f; vertices[2].w = 0.1f;
    vertices[3].x = 0.0f; vertices[3].y = 1.0f; vertices[3].z = 0.0f; vertices[3].w = 0.1f;
    rtcUnmapBuffer(scene,0,RTC_VERTEX_BUFFER);
    AssertNoError();

    int* indices = (int*) rtcMapBuffer(scene,0,RTC_INDEX_BUFFER);
    for (size_t i=0; i<N; i++) {
      indices[i] = 0;
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
    Vertex3fa*   vertices  = (Vertex3fa*  ) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
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
      if (ok0) printf(GREEN("+")); else printf(RED("-"));
      passed &= ok0;
    }
    printf(" %s\n",passed ? GREEN("[PASSED]") : RED("[FAILED]"));
    fflush(stdout);
    numFailedTests += !passed;
  }

  bool rtcore_new_delete_geometry()
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    AssertNoError();
    int geom[128];
    for (size_t i=0; i<128; i++) geom[i] = -1;
    Sphere spheres[128];
    memset(spheres,0,sizeof(spheres));

    for (size_t i=0; i<50; i++) {
      for (size_t j=0; j<10; j++) {
        int index = rand()%128;
        Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
        if (geom[index] == -1) {
          switch (rand()%4) {
          case 0: geom[index] = addSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,10); break;
          case 1: geom[index] = addHair  (scene,RTC_GEOMETRY_STATIC,pos,1.0f,2.0f,10); break;
	  case 2: geom[index] = addSubdivSphere(scene,RTC_GEOMETRY_STATIC,pos,2.0f,4,4); break;
          case 3: 
	    spheres[index] = Sphere(pos,2.0f);
	    geom[index] = addUserGeometryEmpty(scene,&spheres[index]); break;
          }
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
    clearBuffers();
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

  struct RegressionTask
  {
    RegressionTask (size_t sceneIndex, size_t sceneCount, size_t threadCount)
      : sceneIndex(sceneIndex), sceneCount(sceneCount), scene(NULL), numActiveThreads(0) { barrier.init(threadCount); }

    size_t sceneIndex;
    size_t sceneCount;
    RTCScene scene;
    BarrierSys barrier;
    volatile size_t numActiveThreads;
  };

  struct ThreadRegressionTask
  {
    ThreadRegressionTask (size_t threadIndex, size_t threadCount, RegressionTask* task)
      : threadIndex(threadIndex), threadCount(threadCount), task(task) {}

    size_t threadIndex;
    size_t threadCount;
    RegressionTask* task;
  };

  void rtcore_regression_static_thread(void* ptr)
  {
    ThreadRegressionTask* thread = (ThreadRegressionTask*) ptr;
    RegressionTask* task = thread->task;
    if (thread->threadIndex > 0) 
    {
      for (size_t i=0; i<task->sceneCount; i++) 
      {
	task->barrier.wait();
	if (thread->threadIndex < task->numActiveThreads) 
	{
	  rtcCommitThread(task->scene,thread->threadIndex,task->numActiveThreads);
	  CountErrors();
	  
	  for (size_t i=0; i<100; i++)
	    shootRays(task->scene);
	}
	task->barrier.wait();
      }
      delete thread; thread = NULL;
      return;
    }

    CountErrors();
    int geom[1024];
    int types[1024];
    Sphere spheres[1024];
    size_t numVertices[1024];
    for (size_t i=0; i<1024; i++)  {
      geom[i] = -1;
      types[i] = 0;
      numVertices[i] = 0;
    }

    for (size_t i=0; i<task->sceneCount; i++) 
    {
      srand(task->sceneIndex*23565+i*3242);
      if (i%20 == 0) std::cout << "." << std::flush;

      RTCSceneFlags sflag = getSceneFlag(i); 
      task->scene = rtcNewScene(sflag,aflags);
      vector_t<Sphere*> spheres;
      CountErrors();

      for (size_t j=0; j<10; j++) 
      {
        Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
	int type = rand()%6;
#if !defined(__MIC__) 
        switch (rand()%16) {
        case 0: pos = Vec3fa(nan); break;
        case 1: pos = Vec3fa(inf); break;
        case 2: pos = Vec3fa(1E30f); break;
        default: break;
        };
#endif
	size_t numPhi = rand()%100;
	if (type == 2) numPhi = rand()%10;
        size_t numTriangles = 2*2*numPhi*(numPhi-1);
	numTriangles = rand()%(numTriangles+1);
        switch (type) {
        case 0: addSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,0.0f); break;
	case 1: addSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
	case 2: addSubdivSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,4,numTriangles,0.0f); break;
	case 3: addHair  (task->scene,RTC_GEOMETRY_STATIC,pos,1.0f,2.0f,numTriangles,0.0f); break;
	case 4: addHair  (task->scene,RTC_GEOMETRY_STATIC,pos,1.0f,2.0f,numTriangles,1.0f); break; 
        case 5: {
	  Sphere* sphere = new Sphere(pos,2.0f); spheres.push_back(sphere); 
	  addUserGeometryEmpty(task->scene,sphere); break;
        }
	}
        CountErrors();
      }

      if (thread->threadCount) {
	task->numActiveThreads = max(size_t(1),rand() % thread->threadCount);
	task->barrier.wait();
	rtcCommitThread(task->scene,thread->threadIndex,task->numActiveThreads);
      } else {
	rtcCommit(task->scene);
      }
      CountErrors();

      for (size_t i=0; i<100; i++)
        shootRays(task->scene);

      if (thread->threadCount) 
	task->barrier.wait();

      rtcDeleteScene (task->scene);
      CountErrors();

      for (size_t i=0; i<spheres.size(); i++)
	delete spheres[i];
    }
    return;
  }

  void rtcore_regression_dynamic_thread(void* ptr)
  {
    ThreadRegressionTask* thread = (ThreadRegressionTask*) ptr;
    RegressionTask* task = thread->task;
    if (thread->threadIndex > 0) 
    {
      for (size_t i=0; i<task->sceneCount; i++) 
      {
	task->barrier.wait();
	if (thread->threadIndex < task->numActiveThreads) 
	{
	  rtcCommitThread(task->scene,thread->threadIndex,task->numActiveThreads);
	  CountErrors();
	  
	  for (size_t i=0; i<100; i++)
	    shootRays(task->scene);
	}
	task->barrier.wait();
      }
      delete thread; thread = NULL;
      return;
    }
    task->scene = rtcNewScene(RTC_SCENE_DYNAMIC,aflags);
    CountErrors();
    int geom[1024];
    int types[1024];
    Sphere spheres[1024];
    size_t numVertices[1024];
    for (size_t i=0; i<1024; i++)  {
      geom[i] = -1;
      types[i] = 0;
      numVertices[i] = 0;
    }

    for (size_t i=0; i<task->sceneCount; i++) 
    {
      srand(task->sceneIndex*23565+i*3242);
      if (i%20 == 0) std::cout << "." << std::flush;

      for (size_t j=0; j<40; j++) 
      {
        int index = rand()%1024;
        if (geom[index] == -1) 
        {
          int type = rand()%10;
          Vec3fa pos = 100.0f*Vec3fa(drand48(),drand48(),drand48());
#if !defined(__MIC__)
          switch (rand()%16) {
          case 0: pos = Vec3fa(nan); break;
          case 1: pos = Vec3fa(inf); break;
          case 2: pos = Vec3fa(1E30f); break;
          default: break;
          };
#endif
          size_t numPhi = rand()%100;
	  if (type >= 3 || type <= 5) numPhi = rand()%10;
          size_t numTriangles = 2*2*numPhi*(numPhi-1);
          numTriangles = rand()%(numTriangles+1);
          types[index] = type;
          numVertices[index] = 2*numPhi*(numPhi+1);
          switch (type) {
          case 0: geom[index] = addSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,0.0f); break;
          case 1: geom[index] = addSphere(task->scene,RTC_GEOMETRY_DEFORMABLE,pos,2.0f,numPhi,numTriangles,0.0f); break;
          case 2: geom[index] = addSphere(task->scene,RTC_GEOMETRY_DYNAMIC,pos,2.0f,numPhi,numTriangles,0.0f); break;

          case 3: geom[index] = addSubdivSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,4,numTriangles,0.0f); break;
	  case 4: geom[index] = addSubdivSphere(task->scene,RTC_GEOMETRY_DEFORMABLE,pos,2.0f,numPhi,4,numTriangles,0.0f); break;
	  case 5: geom[index] = addSubdivSphere(task->scene,RTC_GEOMETRY_DYNAMIC,pos,2.0f,numPhi,4,numTriangles,0.0f); break;

          case 6: geom[index] = addSphere(task->scene,RTC_GEOMETRY_STATIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
          case 7: geom[index] = addSphere(task->scene,RTC_GEOMETRY_DEFORMABLE,pos,2.0f,numPhi,numTriangles,1.0f); break;
          case 8: geom[index] = addSphere(task->scene,RTC_GEOMETRY_DYNAMIC,pos,2.0f,numPhi,numTriangles,1.0f); break;
            
          case 9: spheres[index] = Sphere(pos,2.0f); geom[index] = addUserGeometryEmpty(task->scene,&spheres[index]); break;
          }; 
	  CountErrors();
        }
        else 
        {
          switch (types[index]) {
          case 0:
          case 3:
          case 6:
	  case 9: {
            rtcDeleteGeometry(task->scene,geom[index]);     
	    CountErrors();
            geom[index] = -1; 
            break;
          }
          case 1: 
          case 2:
          case 4: 
          case 5:
	  case 7: 
          case 8: {
            int op = rand()%2;
            switch (op) {
            case 0: {
              rtcDeleteGeometry(task->scene,geom[index]);     
	      CountErrors();
              geom[index] = -1; 
              break;
            }
            case 1: {
              Vertex3f* vertices = (Vertex3f*) rtcMapBuffer(task->scene,geom[index],RTC_VERTEX_BUFFER);
              for (size_t i=0; i<numVertices[index]; i++) vertices[i] += Vertex3f(0.1f);
              rtcUnmapBuffer(task->scene,geom[index],RTC_VERTEX_BUFFER);

              if (types[index] == 7 || types[index] == 8) {
                Vertex3f* vertices = (Vertex3f*) rtcMapBuffer(task->scene,geom[index],RTC_VERTEX_BUFFER1);
                for (size_t i=0; i<numVertices[index]; i++) vertices[i] += Vertex3f(0.1f);
                rtcUnmapBuffer(task->scene,geom[index],RTC_VERTEX_BUFFER1);
              }
              break;
            }
            }
            break;
          }
          }
        }
      }

      if (thread->threadCount) {
	task->numActiveThreads = max(size_t(1),rand() % thread->threadCount);
	task->barrier.wait();
	rtcCommitThread(task->scene,thread->threadIndex,task->numActiveThreads);
      } else {
	rtcCommit(task->scene);
      }
      CountErrors();

      for (size_t i=0; i<100; i++)
        shootRays(task->scene);

      if (thread->threadCount) 
	task->barrier.wait();
    }

    rtcDeleteScene (task->scene);
    CountErrors();

    delete thread; thread = NULL;
    return;
  }

  bool rtcore_regression (thread_func func, bool userThreads)
  {
    errorCounter = 0;
    size_t sceneIndex = 0;
    while (sceneIndex < regressionN/5) 
    {
      if (userThreads)
      {
	size_t numThreads = getNumberOfLogicalThreads();
#if defined (__MIC__)
	numThreads -= 4;
#endif

        std::vector<RegressionTask*> tasks;
      
	while (numThreads) 
	{
	  size_t N = max(size_t(1),rand()%numThreads); numThreads -= N;
	  RegressionTask* task = new RegressionTask(sceneIndex++,5,N);
	  tasks.push_back(task);

	  for (size_t i=0; i<N; i++) 
	    g_threads.push_back(createThread(func,new ThreadRegressionTask(i,N,task),DEFAULT_STACK_SIZE,numThreads+i));
	}
	
	for (size_t i=0; i<g_threads.size(); i++)
	  join(g_threads[i]);
        for (size_t i=0; i<tasks.size(); i++)
          delete tasks[i];
	
	g_threads.clear();
	clearBuffers();
      }
      else
      {
	func(new ThreadRegressionTask(0,0,new RegressionTask(sceneIndex++,5,0)));
	clearBuffers();
      }	
    }
    return errorCounter == 0;
  }

  bool rtcore_regression_garbage()
  {
    for (size_t i=0; i<5*regressionN; i++) 
    {
      srand(i*23565);
      if (i%20 == 0) std::cout << "." << std::flush;

      RTCSceneFlags sflag = getSceneFlag(i); 
      RTCScene scene = rtcNewScene(sflag,aflags);
      AssertNoError();

      for (size_t j=0; j<20; j++) 
      {
	size_t numTriangles = rand()%256;
        switch (rand()%4) {
        case 0: addGarbageTriangles(scene,RTC_GEOMETRY_STATIC,numTriangles,false); break;
	case 1: addGarbageTriangles(scene,RTC_GEOMETRY_STATIC,numTriangles,true); break;
        case 2: addGarbageHair     (scene,RTC_GEOMETRY_STATIC,numTriangles,false); break;
	case 3: addGarbageHair     (scene,RTC_GEOMETRY_STATIC,numTriangles,true); break;
	}
        AssertNoError();
      }

      rtcCommit(scene);
      AssertNoError();
      rtcDeleteScene (scene);
      clearBuffers();
      AssertNoError();
    }
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
#if 1
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

#if defined(RTCORE_BUFFER_STRIDE)
    POSITIVE("buffer_stride",             rtcore_buffer_stride());
#endif


    POSITIVE("dynamic_enable_disable",    rtcore_dynamic_enable_disable());

    POSITIVE("update_deformable",         rtcore_update(RTC_GEOMETRY_DEFORMABLE));
    POSITIVE("update_dynamic",            rtcore_update(RTC_GEOMETRY_DYNAMIC));
    POSITIVE("overlapping_triangles",     rtcore_overlapping_triangles(100000));
    POSITIVE("overlapping_hair",          rtcore_overlapping_hair(100000));

    rtcore_build();

    POSITIVE("new_delete_geometry",       rtcore_new_delete_geometry());

#if defined(RTCORE_RAY_MASK)
    rtcore_ray_masks_all();
#endif

#if defined(RTCORE_INTERSECTION_FILTER)
    rtcore_filter_all();
#endif

#if defined(RTCORE_BACKFACE_CULLING)
    rtcore_backface_culling_all();
#endif

    rtcore_packet_write_test_all();

    const Vec3fa pos = Vec3fa(148376.0f,1234.0f,-223423.0f);
	rtcore_watertight_closed1("sphere", pos);
    rtcore_watertight_closed1("cube",pos);
    rtcore_watertight_plane1(100000);
#if !defined(__MIC__)
    rtcore_watertight_closed4("sphere",pos);
    rtcore_watertight_closed4("cube",pos);
    rtcore_watertight_plane4(100000);
#endif

#if defined(__TARGET_AVX__) || defined(__TARGET_AVX2__)
    if (has_feature(AVX)) {
      rtcore_watertight_closed8("sphere",pos);
      rtcore_watertight_closed8("cube",pos);
      rtcore_watertight_plane8(100000);
    }
#endif

#if defined(__MIC__)
    rtcore_watertight_closed16("sphere",pos);
    rtcore_watertight_closed16("cube",pos);
    rtcore_watertight_plane16(100000);

#endif

#if defined(__FIX_RAYS__) // FIXME: this mode is deactivated, did never really work
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

    POSITIVE("regression_static",         rtcore_regression(rtcore_regression_static_thread,false));
    POSITIVE("regression_dynamic",        rtcore_regression(rtcore_regression_dynamic_thread,false));

#endif

#if !defined(__MIC__)
    POSITIVE("regression_static_user_threads", rtcore_regression(rtcore_regression_static_thread,true));
    POSITIVE("regression_dynamic_user_threads", rtcore_regression(rtcore_regression_dynamic_thread,true));

    POSITIVE("regression_garbage_geom",   rtcore_regression_garbage());
#endif

    rtcExit();
    return numFailedTests;
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

