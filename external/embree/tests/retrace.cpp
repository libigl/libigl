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
#include "embree2/rtcore.h"
#include "embree2/rtcore_ray.h"
#include "../kernels/common/default.h"
#include "../kernels/common/raystream_log.h"
#include "sys/intrinsics.h"
#include "sys/thread.h"
#include "sys/sysinfo.h"
#include "sys/sync/barrier.h"
#include "sys/sync/mutex.h"
#include "sys/sync/condition.h"
#include <vector>
#include <iostream>
#include <fstream>

#define DBG(x) 

#define SSC_MARK(mark_value)     \
        {__asm  mov ebx, mark_value} \
        {__asm  _emit 0x64}          \
        {__asm  _emit 0x67}          \
        {__asm  _emit 0x90}

namespace embree
{
  struct RayStreamStats
  {
    size_t numTotalRays;
    size_t numRayPackets;
    size_t numIntersectRayPackets;
    size_t numOccludedRayPackets;
    size_t numIntersectRays;
    size_t numOccludedRays;
    size_t num4widePackets;
    size_t num8widePackets;

    RayStreamStats() {
      memset(this,0,sizeof(RayStreamStats));
    }


    void add(RayStreamLogger::LogRay16 &r)
    {
      size_t numRays = r.numRays;
      numRayPackets++;
      numTotalRays += numRays;
      if (r.type == RayStreamLogger::RAY_INTERSECT)
	{
	  numIntersectRayPackets++;
	  numIntersectRays += numRays;
	}
      else if (r.type == RayStreamLogger::RAY_OCCLUDED)
	{
	  numOccludedRayPackets++;
	  numOccludedRays += numRays;
	}
      else
	FATAL("unknown log ray type");
      num4widePackets += (numRays+3)/4;
      num8widePackets += (numRays+7)/8;
    }

    void add(RayStreamLogger::LogRay4 &r)
    {
      size_t numRays = r.numRays;
      numRayPackets++;
      numTotalRays += numRays;
      if (r.type == RayStreamLogger::RAY_INTERSECT)
	{
	  numIntersectRayPackets++;
	  numIntersectRays += numRays;
	}
      else if (r.type == RayStreamLogger::RAY_OCCLUDED)
	{
	  numOccludedRayPackets++;
	  numOccludedRays += numRays;
	}
      else
	FATAL("unknown log ray type");
    }

    void add(RayStreamLogger::LogRay8 &r)
    {
      size_t numRays = r.numRays;
      numRayPackets++;
      numTotalRays += numRays;
      if (r.type == RayStreamLogger::RAY_INTERSECT)
	{
	  numIntersectRayPackets++;
	  numIntersectRays += numRays;
	}
      else if (r.type == RayStreamLogger::RAY_OCCLUDED)
	{
	  numOccludedRayPackets++;
	  numOccludedRays += numRays;
	}
      else
	FATAL("unknown log ray type");
      num4widePackets += (numRays+3)/4;
    }

    void add(RayStreamLogger::LogRay1 &r)
    {
      numRayPackets++;
      numTotalRays ++;
      if (r.type == RayStreamLogger::RAY_INTERSECT)
	{
	  numIntersectRayPackets++;
	  numIntersectRays++;
	}
      else if (r.type == RayStreamLogger::RAY_OCCLUDED)
	{
	  numOccludedRayPackets++;
	  numOccludedRays++;
	}
      else
	FATAL("unknown log ray type");
      num4widePackets += 1;
      num8widePackets += 1;
    }

    void print(size_t simd_width)
    {
      std::cout << "numTotalRays                      = " << numTotalRays << std::endl;
      std::cout << "numRayPackets                     = " << numRayPackets << std::endl;
      std::cout << "numIntersectionRays               = " << numIntersectRays << " [" << 100. * (double)numIntersectRays / numTotalRays << "%]" << std::endl;
      std::cout << "numOcclusionRays                  = " << numOccludedRays << " [" << 100. * (double)numOccludedRays / numTotalRays << "%]" << std::endl;
      if (simd_width > 1)
        {
          std::cout << "avg. intersect " << simd_width << "-wide packet utilization = " << 100. *  (double)numIntersectRays / (numIntersectRayPackets * (double)simd_width) << "%" << std::endl;
          std::cout << "avg. occluded " << simd_width << "-wide packet utilization = " << 100. *  (double)numOccludedRays  / (numOccludedRayPackets  * (double)simd_width) << "%" << std::endl;
          std::cout << "avg. total " << simd_width << "-wide packet utilization     = " << 100. * (double)numTotalRays / (numRayPackets * (double)simd_width)  << "%" << std::endl;
        }
      if (simd_width == 16)
        {
          std::cout << "avg. 4-wide packet utilization    = " << 100. * (double)numTotalRays / (num4widePackets * 4.)  << "%" << std::endl;
          std::cout << "avg. 8-wide packet utilization    = " << 100. * (double)numTotalRays / (num8widePackets * 8.)  << "%" << std::endl;
        }
    } 
  };


  struct RetraceTask
  {
    RTCScene scene;
    void *raydata;
    void *raydata_verify;
    size_t numLogRayStreamElements;
    bool check;
  };  

  /* configuration */

#if !defined(__MIC__)
  RTCAlgorithmFlags aflags = (RTCAlgorithmFlags) (RTC_INTERSECT1 | RTC_INTERSECT4 | RTC_INTERSECT8);
  static std::string g_rtcore = "verbose=2";
#else
  RTCAlgorithmFlags aflags = (RTCAlgorithmFlags) (RTC_INTERSECT1 | RTC_INTERSECT16);
  static std::string g_rtcore = "verbose=2,threads=1";
#endif

  /* vertex and triangle layout */
  struct Vertex   { float x,y,z,a; };
  struct Triangle { int v0, v1, v2; };

  __forceinline std::ostream &operator<<(std::ostream &o, const Vertex &v)
  {
    o << "vtx " << v.x << " " << v.y << " " << v.z << " " << v.a << std::endl;
    return o;
  } 

  __forceinline std::ostream &operator<<(std::ostream &o, const Triangle &t)
  {
    o << "tri " << t.v0 << " " << t.v1 << " " << t.v2 << std::endl;
    return o;
  } 

  static AlignedAtomicCounter32 g_counter = 0;
  static bool g_check = false;
  static bool g_sde = false;
  static size_t g_threadCount = 1;
  static size_t g_frames = 1;
  static size_t g_simd_width = 0;
  static AlignedAtomicCounter32 g_rays_traced = 0;
  static AlignedAtomicCounter32 g_rays_traced_diff = 0;
  static std::vector<thread_t> g_threads;
#if !defined(__MIC__)
  static BarrierSys g_barrier;
#else
  static LinearBarrierActive g_barrier;
#endif

  static bool g_exitThreads = false;
  static RetraceTask g_retraceTask;
  static MutexSys g_mutex;

    
#if defined(__MIC__)
  static std::string g_binaries_path = "/home/micuser/";
#else
  static std::string g_binaries_path = "./";
#endif


#define AssertNoError()                                 \
  if (rtcGetError() != RTC_NO_ERROR) return false;
#define AssertAnyError()                                \
  if (rtcGetError() == RTC_NO_ERROR) return false;
#define AssertError(code)                       \
  if (rtcGetError() != code) return false;

  static void parseCommandLine(int argc, char** argv)
  {
    for (int i=1; i<argc; i++)
      {
        std::string tag = argv[i];
        if (tag == "") return;

        else if (tag == "-rtcore" && i+1<argc) {
          g_rtcore += std::stringOf(',') + argv[++i];
        }
        /* rtcore configuration */
        else if (tag == "-check") {
          g_check = true;
        }      
        else if (tag == "-threads" && i+1<argc) {
          g_threadCount = atoi(argv[++i]);
        }
        else if (tag == "-frames" && i+1<argc) {
          g_frames = atoi(argv[++i]);
        }
        else if (tag == "-simd_width" && i+1<argc) {
          g_simd_width = atoi(argv[++i]);
          if (g_simd_width != 1 && g_simd_width != 4 && g_simd_width != 8 && g_simd_width != 16)
            std::cout << "only simd widths of 1,4,8, and 16 are supported" << std::endl;
        }
        else if (tag == "-sde") {
          g_sde = true;
        }
        else if (tag == "-h" || tag == "-help") {
          std::cout << "Usage: retrace [OPTIONS] [PATH_TO_BINARY_FILES] " << std::endl;
          std::cout << "Options:" << std::endl;
          std::cout << "-threads N   : sets number of render/worker threads for the retracing phase to N" << std::endl;
          std::cout << "-frames N    : retraces all rays N times  " << std::endl;
          std::cout << "-check       : loads second ray stream file and validates result of rtcIntersectN/rtcOccludedN for each ray/packet" << std::endl;
          std::cout << "-simd_width N: loads ray stream for simd width N (if existing)" << std:: endl;
          std::cout << "-sde         : inserts markers for generating instruction traces with SDE" << std:: endl;
          exit(0);
        }

        /* skip unknown command line parameter */
        else {

          g_binaries_path = tag;
        }
      }
  }


  bool existsFile(std::string &filename)
  {
    std::ifstream file;
    file.open(filename.c_str(),std::ios::in | std::ios::binary);
    if (!file) return false;
    file.close();
    return true;    
  }

  void *loadGeometryData(std::string &geometryFile)
  {
    std::ifstream geometryData;

    geometryData.open(geometryFile.c_str(),std::ios::in | std::ios::binary);
    if (!geometryData) { FATAL("could not open geometry data file"); }
    geometryData.seekg(0, std::ios::beg);
    std::streampos begin = geometryData.tellg();
    geometryData.seekg(0, std::ios::end);
    std::streampos end = geometryData.tellg();
    size_t fileSize = end - begin;
    char *ptr = (char*)os_malloc(fileSize);
    geometryData.seekg(0, std::ios::beg);
    geometryData.read(ptr,fileSize);
    geometryData.close();
    return ptr;
  }

  template<class T>
  void *loadRayStreamData(std::string &rayStreamFile, size_t &numLogRayStreamElements)
  {
    std::ifstream rayStreamData;

    rayStreamData.open(rayStreamFile.c_str(),std::ios::in | std::ios::binary);
    if (!rayStreamData) { FATAL("could not open raystream data file"); }
    rayStreamData.seekg(0, std::ios::beg);
    std::streampos begin = rayStreamData.tellg();
    rayStreamData.seekg(0, std::ios::end);
    std::streampos end = rayStreamData.tellg();
    size_t fileSize = end - begin;
    char *ptr = (char*)os_malloc(fileSize);
    rayStreamData.seekg(0, std::ios::beg);
    rayStreamData.read(ptr,fileSize);
    numLogRayStreamElements = fileSize / sizeof(T);
    rayStreamData.close();
    return ptr;
  }

  template<class T>
  RayStreamStats analyseRayStreamData(void *r, size_t numLogRayStreamElements)
  {
    RayStreamStats stats;
    std::cout << "numLogRayStreamElements " << numLogRayStreamElements << std::endl;
    for (size_t i=0;i<numLogRayStreamElements;i++)
      stats.add(((T*)r)[i]);
    return stats;
  }

  RTCScene transferGeometryData(char *g)
  {
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,aflags);

    int magick = *(size_t*)g; g += sizeof(int);
    if (magick != 0x35238765LL) {
      FATAL("invalid binary file");
    }

    int numGroups = *(int*)g; g += sizeof(int);

    for (size_t i=0; i<numGroups; i++) 
    {
      int type = *(int*)g; g += sizeof(int);

      if (type == 1)
      {
	int numTimeSteps = *(int*)g; g += sizeof(int);
	int numVertices  = *(int*)g; g += sizeof(int);
	int numTriangles = *(int*)g; g += sizeof(int);
	unsigned int geometry = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, numTriangles, numVertices, numTimeSteps);

	for (size_t i=0; i<numTimeSteps; i++) {
	  if (((size_t)g % 16) != 0) g += 16 - ((size_t)g % 16);
	  rtcSetBuffer(scene, geometry, (RTCBufferType)(RTC_VERTEX_BUFFER0+i), g, 0, sizeof(Vec3fa));
	  g += numVertices*sizeof(Vec3fa);
	}

	if (((size_t)g % 16) != 0) g += 16 - ((size_t)g % 16);
	rtcSetBuffer(scene, geometry, RTC_INDEX_BUFFER,  g, 0, sizeof(Triangle));
	g += numTriangles*sizeof(Triangle);
      }

      else if (type == 2)
      {
	int numTimeSteps = *(int*)g; g += sizeof(int);
	int numVertices  = *(int*)g; g += sizeof(int);
	int numCurves    = *(int*)g; g += sizeof(int);

	unsigned int geometry = rtcNewHairGeometry (scene, RTC_GEOMETRY_STATIC, numCurves, numVertices, numTimeSteps);
	
	for (size_t i=0; i<numTimeSteps; i++) {
	  if (((size_t)g % 16) != 0) g += 16 - ((size_t)g % 16);
	  rtcSetBuffer(scene, geometry, (RTCBufferType)(RTC_VERTEX_BUFFER0+i), g, 0, sizeof(Vec3fa));
	  g += numVertices*sizeof(Vec3fa);
	}

	if (((size_t)g % 16) != 0) g += 16 - ((size_t)g % 16);
	rtcSetBuffer(scene, geometry, RTC_INDEX_BUFFER,  g, 0, sizeof(int));
	g += numCurves*sizeof(int);
      }

      else if (type == -1) {
      }

      else {
	FATAL("unknown geometry type");
      }
    }

    rtcCommit(scene);
    return scene;
  }

  size_t check_ray1_packets(RTCRay &start, RTCRay &end)
  {
    if (start.primID != end.primID) return 1;
    if (start.geomID != end.geomID) return 1;
    if (start.u != end.u) return 1;
    if (start.v != end.v) return 1;
    if (start.tfar != end.tfar) return 1;
    return 0;
  }

  template<class T>
  void print(const char *name, T *ptr, const size_t N)
  {
    std::cout << name << " "; for (size_t i=0;i<N;i++) std::cout << ptr[i] << " "; std::cout << std::endl;
  }

  template<class T>  
  void print_packet(T &t)
  {
    const size_t elements = sizeof(T) / (18 * 4); // 18 elements, and 4 bytes per element
    print("orgx"  ,t.orgx,elements);
    print("orgy"  ,t.orgy,elements);
    print("orgz"  ,t.orgz,elements);
    print("dirx"  ,t.dirx,elements);
    print("diry"  ,t.diry,elements);
    print("dirz"  ,t.dirz,elements);
    print("tnear" ,t.tnear,elements);
    print("tfar"  ,t.tfar,elements);
    print("primID",t.primID,elements);
    print("geomID",t.geomID,elements);
    print("u"     ,t.u,elements);
    print("v"     ,t.v,elements);
    print("Ngx"   ,t.Ngx,elements);
    print("Ngy"   ,t.Ngy,elements);
    print("Ngz"   ,t.Ngz,elements);    
  }

  template<class T>
  size_t check_ray_packets(const unsigned int m_valid, T &start, T &end)
  {
    size_t diff = 0;
    const size_t elements = sizeof(T) / (18 * 4); // 18 elements, and 4 bytes per element

    for (size_t i=0;i<elements;i++)
      {
        if ( (((unsigned int)1 << i) & m_valid) == 0) continue;
        if (start.primID[i] != end.primID[i]) { diff++; continue; }
        if (start.geomID[i] != end.geomID[i]) { diff++; continue; }
        if (start.u[i]      != end.u[i])      { diff++; continue; }
        if (start.v[i]      != end.v[i])      { diff++; continue; }
        if (start.tnear[i]  != end.tnear[i])  { diff++; continue; }
        if (start.tfar[i]   != end.tfar[i])   { diff++; continue; }
        if (start.Ngx[i]    != end.Ngx[i])    { diff++; continue; }
        if (start.Ngy[i]    != end.Ngy[i])    { diff++; continue; }
        if (start.Ngz[i]    != end.Ngz[i])    { diff++; continue; }
      }
    return diff;
  }


#define RAY_BLOCK_SIZE 16

  template<size_t SIMD_WIDTH>
  void retrace_loop()
  {
    size_t rays = 0;
    size_t diff = 0;
    
    while(1)
      {
	size_t global_index = g_counter.add(RAY_BLOCK_SIZE);
	if (global_index >= g_retraceTask.numLogRayStreamElements) break;
	size_t startID = global_index;
	size_t endID   = min(g_retraceTask.numLogRayStreamElements,startID+RAY_BLOCK_SIZE);

	for (size_t index=startID;index<endID;index++)
	  {
            if (SIMD_WIDTH == 1)
              {
                RayStreamLogger::LogRay1 *raydata        = (RayStreamLogger::LogRay1 *)g_retraceTask.raydata;
                RayStreamLogger::LogRay1 *raydata_verify = (RayStreamLogger::LogRay1 *)g_retraceTask.raydata_verify;

                RTCRay &ray = raydata[index].ray;
                rays ++;
                if (raydata[index].type == RayStreamLogger::RAY_INTERSECT)
                  rtcIntersect(g_retraceTask.scene,ray);
                else 
                  rtcOccluded(g_retraceTask.scene,ray);

                if (unlikely(g_check))
                  diff += check_ray1_packets(ray, raydata_verify[index].ray);                    
              }
#if !defined(__MIC__)
            else if (SIMD_WIDTH == 4)
              {
                RayStreamLogger::LogRay4 *raydata        = (RayStreamLogger::LogRay4 *)g_retraceTask.raydata;
                RayStreamLogger::LogRay4 *raydata_verify = (RayStreamLogger::LogRay4 *)g_retraceTask.raydata_verify;

                RTCRay4 &ray4 = raydata[index].ray4;
                sseb valid((int)raydata[index].m_valid);
                rays += raydata[index].numRays;
                if (raydata[index].type == RayStreamLogger::RAY_INTERSECT)
                  rtcIntersect4(&valid,g_retraceTask.scene,ray4);
                else 
                  rtcOccluded4(&valid,g_retraceTask.scene,ray4);

                if (unlikely(g_check))
                  diff += check_ray_packets<RTCRay4>(raydata[index].m_valid,  raydata[index].ray4, raydata_verify[index].ray4);
              }
            else if (SIMD_WIDTH == 8)
              {
                RayStreamLogger::LogRay8 *raydata        = (RayStreamLogger::LogRay8 *)g_retraceTask.raydata;
                RayStreamLogger::LogRay8 *raydata_verify = (RayStreamLogger::LogRay8 *)g_retraceTask.raydata_verify;

                RTCRay8 &ray8 = raydata[index].ray8;
                __aligned(64) sseb valid[2];

                valid[0] = sseb((int)(raydata[index].m_valid & 0xf));
                valid[1] = sseb((int)(raydata[index].m_valid>>4));

                rays += raydata[index].numRays;

                if (raydata[index].type == RayStreamLogger::RAY_INTERSECT)
                  rtcIntersect8(valid,g_retraceTask.scene,ray8);
                else 
                  rtcOccluded8(valid,g_retraceTask.scene,ray8);

                if (unlikely(g_check))
                  diff += check_ray_packets<RTCRay8>(raydata[index].m_valid,  raydata[index].ray8, raydata_verify[index].ray8);
              }
#endif
            else if (SIMD_WIDTH == 16)
              {
                RayStreamLogger::LogRay16 *raydata        = (RayStreamLogger::LogRay16 *)g_retraceTask.raydata;
                RayStreamLogger::LogRay16 *raydata_verify = (RayStreamLogger::LogRay16 *)g_retraceTask.raydata_verify;
#if defined(__MIC__)
                RTCRay16 &ray16 = raydata[index].ray16;
                mic_i valid = select((mic_m)raydata[index].m_valid,mic_i(-1),mic_i(0));
                rays += raydata[index].numRays;

                //raydata[index+1].prefetchL2();
  
                if (raydata[index].type == RayStreamLogger::RAY_INTERSECT)
                  rtcIntersect16(&valid,g_retraceTask.scene,ray16);
                else 
                  rtcOccluded16(&valid,g_retraceTask.scene,ray16);
#endif
                if (unlikely(g_check))
                  diff += check_ray_packets<RTCRay16>(raydata[index].m_valid,  raydata[index].ray16, raydata_verify[index].ray16);
              }


	  }
      }
    if (unlikely(g_check && diff))
      g_rays_traced_diff.add(diff);
    g_rays_traced.add(rays);
  }



  void renderMainLoop(size_t id)
  {
    switch(g_simd_width)
      {
      case 1:
        retrace_loop<1>();
        break;
      case 4:
        retrace_loop<4>();
        break;
      case 8:
        retrace_loop<8>();
        break;
      case 16:
        retrace_loop<16>();
        break;
      };
  }

  void threadMainLoop(void *ptr)
  {
    size_t id = (size_t)ptr;
    setAffinity(id);

    DBG(
        g_mutex.lock();
        DBG_PRINT(id);
        g_mutex.unlock();
        );

    while(1)
      {
        g_barrier.wait(id,g_threadCount);
    
        if (g_exitThreads) break;

        renderMainLoop(id);

        g_barrier.wait(id,g_threadCount);
      }
  }

  void createThreads(size_t numThreads)
  {
    for (size_t i=1; i<numThreads; i++)
      g_threads.push_back(createThread(threadMainLoop,(void*)i,1000000,i));
  }

  /* main function in embree namespace */
  int main(int argc, char** argv) 
  {
    g_threadCount = getNumberOfLogicalThreads(); 
#if defined (__MIC__)
    g_threadCount -= 4;
#endif

    setAffinity(0);
    /* parse command line */  
    parseCommandLine(argc,argv);

    /* perform tests */
    DBG_PRINT(g_rtcore.c_str());
    rtcInit(g_rtcore.c_str());

    DBG_PRINT(g_threadCount);

    /* binary file path */
    g_binaries_path += "/";
    std::cout << "binary file path = " <<  g_binaries_path << std::endl;

    /* load geometry file */
    std::string geometryFileName = g_binaries_path + "geometry.bin";
    std::cout << "loading geometry data from file '" << geometryFileName << "'..." << std::flush;    
    void *g = loadGeometryData(geometryFileName);
    std::cout <<  "done" << std::endl << std::flush;

    /* transfer geometry data */
    std::cout << "transfering geometry data:" << std::endl << std::flush;
    RTCScene scene = transferGeometryData((char*)g);

    /* looking for ray stream file */
    std::string rayStreamFileName;
    std::string rayStreamVerifyFileName;

    if (g_simd_width != 0)
      {
        rayStreamFileName = g_binaries_path + "ray" + std::stringOf(g_simd_width) + ".bin";
        rayStreamVerifyFileName = g_binaries_path + "ray" + std::stringOf(g_simd_width) + "_verify.bin";
      }
    else
      {
        /* looking for stream files in the following order: ray1.bin, ray4.bin, ray8.bin, ray16.bin */
        for (size_t shift=0;shift<=4;shift++)
          {
            g_simd_width = (size_t)1 << shift;
            rayStreamFileName = g_binaries_path + "ray" + std::stringOf(g_simd_width) + ".bin";
            rayStreamVerifyFileName = g_binaries_path + "ray" + std::stringOf(g_simd_width) + "_verify.bin";
            if (existsFile( rayStreamFileName )) break;

          }
      }
   
    if (g_simd_width == 0)
      FATAL("no valid ray stream data files found");

    DBG_PRINT( rayStreamFileName );
    DBG_PRINT( rayStreamVerifyFileName );

    if (!existsFile( rayStreamFileName )) FATAL("ray stream file does not exists!");
    if (!existsFile( rayStreamVerifyFileName )) FATAL("ray stream verify file does not exists!");


    /* load ray stream data */
    std::cout << "loading ray stream data from files '" << rayStreamFileName << "/" << rayStreamVerifyFileName << "'..." << std::flush;    
    size_t numLogRayStreamElements       = 0;
    size_t numLogRayStreamElementsVerify = 0;

    void *raydata        = NULL;
    void *raydata_verify = NULL;

    switch(g_simd_width)
      {
      case 1:
        raydata = loadRayStreamData<RayStreamLogger::LogRay1>(rayStreamFileName, numLogRayStreamElements);
        if (g_check)
          raydata_verify = loadRayStreamData<RayStreamLogger::LogRay1>(rayStreamVerifyFileName, numLogRayStreamElementsVerify); 
        break;
      case 4:
        raydata = loadRayStreamData<RayStreamLogger::LogRay4>(rayStreamFileName, numLogRayStreamElements);
        if (g_check)
          raydata_verify = loadRayStreamData<RayStreamLogger::LogRay4>(rayStreamVerifyFileName, numLogRayStreamElementsVerify); 
        break;
      case 8:
        raydata = loadRayStreamData<RayStreamLogger::LogRay8>(rayStreamFileName, numLogRayStreamElements);
        if (g_check)
          raydata_verify = loadRayStreamData<RayStreamLogger::LogRay8>(rayStreamVerifyFileName, numLogRayStreamElementsVerify); 
        break;
      case 16:
        raydata = loadRayStreamData<RayStreamLogger::LogRay16>(rayStreamFileName, numLogRayStreamElements);
        if (g_check)
          raydata_verify = loadRayStreamData<RayStreamLogger::LogRay16>(rayStreamVerifyFileName, numLogRayStreamElementsVerify); 
        break;
      default:
        FATAL("unknown SIMD width");
      }

    std::cout <<  "done" << std::endl << std::flush;

    if (g_check)
      if (numLogRayStreamElements != numLogRayStreamElementsVerify)
        FATAL("numLogRayStreamElements != numLogRayStreamElementsVerify");

    /* analyse ray stream data */
    std::cout << "analyse ray stream:" << std::endl << std::flush;    
    RayStreamStats stats;
    switch(g_simd_width)
      {
      case 1:
        stats = analyseRayStreamData<RayStreamLogger::LogRay1>(raydata,numLogRayStreamElements);
        break;
      case 4:
        stats = analyseRayStreamData<RayStreamLogger::LogRay4>(raydata,numLogRayStreamElements);
        break;
      case 8:
        stats = analyseRayStreamData<RayStreamLogger::LogRay8>(raydata,numLogRayStreamElements);
        break;
      case 16:
        stats = analyseRayStreamData<RayStreamLogger::LogRay16>(raydata,numLogRayStreamElements);
        break;
      }

    stats.print(g_simd_width);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    FATAL("ray stream logger still active, must be disabled to run 'retrace'");
#endif



    /* init global tasking barrier */
    g_barrier.init( g_threadCount );


    g_retraceTask.scene                   = scene;
    g_retraceTask.raydata                 = raydata;
    g_retraceTask.raydata_verify          = raydata_verify;
    g_retraceTask.numLogRayStreamElements = numLogRayStreamElements;
    g_retraceTask.check                   = g_check;


    std::cout << "using " << g_threadCount << " threads for retracing rays" << std::endl << std::flush;
    createThreads(g_threadCount);


    /* retrace ray packets */
    DBG_PRINT( g_threadCount );

    std::cout << "Retracing logged rays:" << std::endl << std::flush;
    
    double avg_time = 0;
    double mrays_sec = 0;
    for (size_t i=0;i<g_frames;i++)
      {
	double dt = getSeconds();

        /* retrace rays using all threads */

	g_rays_traced = 0;
	g_rays_traced_diff = 0;
	g_counter = 0;

	if (g_sde)
	  {
#if defined(__INTEL_COMPILER)
	    __asm { int 3 };
	    SSC_MARK(111);
#endif
	  }

        g_barrier.wait(0,g_threadCount);
        renderMainLoop(0);

	dt = getSeconds()-dt;
	mrays_sec += (double)g_rays_traced / dt / 1000000.;
#if 0
	std::cout << "frame " << i << " => time " << 1000. * dt << " " << 1. / dt << " fps " << "ms " << g_rays_traced / dt / 1000000. << " mrays/sec" << std::endl;
#endif

        g_barrier.wait(0,g_threadCount);

	if (g_sde)
	  {
#if defined(__INTEL_COMPILER)
	    __asm { int 3 };
	    SSC_MARK(222);
#endif
	  }

        if (unlikely(g_check))
          std::cout << g_rays_traced_diff << " rays differ in result (" << 100. * g_rays_traced_diff / g_rays_traced << "%)" << std::endl;
      }
    std::cout << "rays " << g_rays_traced << " avg. mrays/sec = " << mrays_sec / (double)g_frames << std::endl;

    std::cout << "freeing threads..." << std::flush;
    g_exitThreads = true;
    g_barrier.wait(0,g_threadCount);
    std::cout << "done" << std::endl << std::flush;


    /* done */
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
