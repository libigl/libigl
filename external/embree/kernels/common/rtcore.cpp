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

#ifdef _WIN32
#  define RTCORE_API extern "C" __declspec(dllexport)
#else
#  define RTCORE_API extern "C" __attribute__ ((visibility ("default")))
#endif

#include "common/default.h"
#include "common/alloc.h"
#include "embree2/rtcore.h"
#include "common/scene.h"
#include "sys/taskscheduler.h"
#include "sys/thread.h"
#include "raystream_log.h"

#define TRACE(x) //std::cout << #x << std::endl;

namespace embree
{
#define CATCH_BEGIN try {
#define CATCH_END                                                       \
  } catch (std::bad_alloc&) {                                           \
    process_error(RTC_OUT_OF_MEMORY,"out of memory");                   \
  } catch (std::exception& e) {                                         \
    process_error(RTC_UNKNOWN_ERROR,e.what());                          \
 } catch (...) {                                                        \
    process_error(RTC_UNKNOWN_ERROR,"unknown exception caught");        \
  }

#define VERIFY_HANDLE(handle) \
  if (handle == NULL) {                                                 \
    process_error(RTC_INVALID_ARGUMENT,"invalid argument");             \
  }

#define VERIFY_GEOMID(id) \
  if (id == -1) {                                                 \
    process_error(RTC_INVALID_ARGUMENT,"invalid argument");       \
  }
  
  /* functions to initialize global state */
  void init_globals();

  /* register functions for accels */
  void BVH4Register();
  void BVH8Register();

#if defined(__MIC__)
  void BVH4iRegister();
  void BVH4MBRegister();
  void BVH4HairRegister();
#endif

  /*! intersector registration functions */
  DECLARE_SYMBOL(RTCBoundsFunc,InstanceBoundsFunc);
  DECLARE_SYMBOL(AccelSet::Intersector1,InstanceIntersector1);
  DECLARE_SYMBOL(AccelSet::Intersector4,InstanceIntersector4);
  DECLARE_SYMBOL(AccelSet::Intersector8,InstanceIntersector8);
  DECLARE_SYMBOL(AccelSet::Intersector16,InstanceIntersector16);
  
  /* global settings */
  std::string g_tri_accel = "default";                 //!< acceleration structure to use for triangles
  std::string g_tri_builder = "default";               //!< builder to use for triangles
  std::string g_tri_traverser = "default";             //!< traverser to use for triangles
  double      g_tri_builder_replication_factor = 2.0f; //!< maximally factor*N many primitives in accel

  std::string g_tri_accel_mb = "default";              //!< acceleration structure to use for motion blur triangles
  std::string g_tri_builder_mb = "default";            //!< builder to use for motion blur triangles
  std::string g_tri_traverser_mb = "default";          //!< traverser to use for triangles

  std::string g_hair_accel = "default";                //!< hair acceleration structure to use
  std::string g_hair_builder = "default";              //!< builder to use for hair
  std::string g_hair_traverser = "default";             //!< traverser to use for hair
  double      g_hair_builder_replication_factor = 3.0f; //!< maximally factor*N many primitives in accel
  float       g_memory_preallocation_factor = 1.0f; 

  std::string g_subdiv_accel = "default";               //!< acceleration structure to use for subdivision surfaces

  int g_scene_flags = -1;                               //!< scene flags to use
  size_t g_verbose = 0;                                 //!< verbosity of output
  size_t g_numThreads = 0;                              //!< number of threads to use in builders
  size_t g_benchmark = 0;
  size_t g_regression_testing = 0;                      //!< enables regression tests at startup

  void initSettings()
  {
    g_tri_accel = "default";
    g_tri_builder = "default";
    g_tri_traverser = "default";
    g_tri_builder_replication_factor = 2.0f;

    g_tri_accel_mb = "default";
    g_tri_builder_mb = "default";
    g_tri_traverser_mb = "default";

    g_hair_accel = "default";
    g_hair_builder = "default";
    g_hair_traverser = "default";
    g_hair_builder_replication_factor = 3.0f;
    g_memory_preallocation_factor = 1.0f;

    g_subdiv_accel = "default";

    g_scene_flags = -1;
    g_verbose = 0;
    g_numThreads = 0;
    g_benchmark = 0;
  }

  void printSettings()
  {
    std::cout << "general:" << std::endl;
    std::cout << "  build threads = " << g_numThreads << std::endl;
    std::cout << "  verbosity     = " << g_verbose << std::endl;

    std::cout << "triangles:" << std::endl;
    std::cout << "  accel         = " << g_tri_accel << std::endl;
    std::cout << "  builder       = " << g_tri_builder << std::endl;
    std::cout << "  traverser     = " << g_tri_traverser << std::endl;
    std::cout << "  replications  = " << g_tri_builder_replication_factor << std::endl;

    std::cout << "motion blur triangles:" << std::endl;
    std::cout << "  accel         = " << g_tri_accel_mb << std::endl;
    std::cout << "  builder       = " << g_tri_builder_mb << std::endl;
    std::cout << "  traverser     = " << g_tri_traverser_mb << std::endl;

    std::cout << "hair:" << std::endl;
    std::cout << "  accel         = " << g_hair_accel << std::endl;
    std::cout << "  builder       = " << g_hair_builder << std::endl;
    std::cout << "  traverser     = " << g_hair_traverser << std::endl;
    std::cout << "  replications  = " << g_hair_builder_replication_factor << std::endl;

    std::cout << "subdivision surfaces:" << std::endl;
    std::cout << "  accel         = " << g_subdiv_accel << std::endl;

#if defined(__MIC__)
    std::cout << "memory allocation:" << std::endl;
    std::cout << "  preallocation_factor  = " << g_memory_preallocation_factor << std::endl;
#endif
  }
  
  /* error flag */
  static tls_t g_error = NULL;
  static std::vector<RTCError*> g_errors;
  static MutexSys g_errors_mutex;
  static RTC_ERROR_FUNCTION g_error_function = NULL;

  /* mutex to make API thread safe */
  static MutexSys g_mutex;

  /* set if embree got initialized */
  static bool g_initialized = false;

  void skipSpace(const char* str, size_t& pos) {
    while (str[pos] == ' ') pos++;
  }

  int parseInt(const char* str, size_t& pos) 
  {
    skipSpace(str,pos);
    size_t begin = pos;
    while (isdigit(str[pos])) pos++;
    return atoi(str+begin);
  }

  float parseFloat(const char* str, size_t& pos) 
  {
    skipSpace(str,pos);
    size_t begin = pos;
    while (isdigit(str[pos]) || str[pos] == '.') pos++;
    return atof(str+begin);
  }

  std::string parseIdentifier(const char* str, size_t& pos) 
  {
    skipSpace(str,pos);
    size_t begin = pos;
    while (isalnum(str[pos]) || str[pos] == '_' || str[pos] == '.') pos++;
    return std::string(str+begin,str+pos);
  }

  bool parseSymbol(const char* str, char c, size_t& pos) 
  {
    skipSpace(str,pos);
    if (str[pos] == c) { pos++; return true; }
    return false;
  }

  bool findNext(const char* str, char c, size_t& pos) 
  {
    while (str[pos] && str[pos] != c) pos++;
    if (str[pos] == c) { pos++; return true; }
    else return false;
  }

  void InstanceIntersectorsRegister ()
  {
    int features = getCPUFeatures();
#if defined(__MIC__)
    SELECT_SYMBOL_KNC(features,InstanceBoundsFunc);
    SELECT_SYMBOL_KNC(features,InstanceIntersector1);
    SELECT_SYMBOL_KNC(features,InstanceIntersector16);
#else
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,InstanceBoundsFunc);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,InstanceIntersector1);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,InstanceIntersector4);
    SELECT_SYMBOL_AVX_AVX2(features,InstanceIntersector8);
#endif
  }

  LockStepTaskScheduler regression_task_scheduler;

  void task_regression_testing(void* This, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* taskGroup) 
  {
    if (regression_tests == NULL) return;
    LockStepTaskScheduler::setInstance(&regression_task_scheduler);
    LockStepTaskScheduler::Init init(threadIndex,threadCount,&regression_task_scheduler);
    if (threadIndex != 0) return;
    for (size_t i=0; i<regression_tests->size(); i++) 
      (*(*regression_tests)[i])();
  }

  RTCORE_API void rtcInit(const char* cfg) 
  {
    Lock<MutexSys> lock(g_mutex);
    TRACE(rtcInit);
    CATCH_BEGIN;

    if (g_initialized) {
      g_mutex.unlock();
      process_error(RTC_INVALID_OPERATION,"already initialized");
      g_mutex.lock();
      return;
    }
    g_initialized = true;

    /* reset global state */
    initSettings();
    
    if (cfg != NULL) 
    {
      size_t pos = 0;
      do {
        std::string tok = parseIdentifier (cfg,pos);

        if (tok == "threads" && parseSymbol(cfg,'=',pos)) 
	{
	  g_numThreads = parseInt(cfg,pos);
#if defined(__MIC__)
	  if (!(g_numThreads == 1 || (g_numThreads % 4) == 0)) {
	    g_mutex.unlock();
	    process_error(RTC_INVALID_OPERATION,"Xeon Phi supports only number of threads % 4 == 0, or threads == 1");
	    g_mutex.lock();
            return;
          }
#endif
        }
        else if (tok == "isa" && parseSymbol (cfg,'=',pos)) 
	{
	  std::string isa = parseIdentifier (cfg,pos);
	  if      (isa == "sse" ) cpu_features = SSE;
	  else if (isa == "sse2") cpu_features = SSE2;
	  else if (isa == "sse3") cpu_features = SSE3;
	  else if (isa == "ssse3") cpu_features = SSSE3;
	  else if (isa == "sse41") cpu_features = SSE41;
	  else if (isa == "sse4.1") cpu_features = SSE41;
	  else if (isa == "sse42") cpu_features = SSE42;
	  else if (isa == "sse4.2") cpu_features = SSE42;
	  else if (isa == "avx") cpu_features = AVX;
	  else if (isa == "avxi") cpu_features = AVXI;
	  else if (isa == "avx2") cpu_features = AVX2;
	}

        else if ((tok == "tri_accel" || tok == "accel") && parseSymbol (cfg,'=',pos))
            g_tri_accel = parseIdentifier (cfg,pos);
	else if ((tok == "tri_builder" || tok == "builder") && parseSymbol (cfg,'=',pos))
	    g_tri_builder = parseIdentifier (cfg,pos);
	else if ((tok == "tri_traverser" || tok == "traverser") && parseSymbol (cfg,'=',pos))
            g_tri_traverser = parseIdentifier (cfg,pos);
	else if (tok == "tri_builder_replication_factor" && parseSymbol (cfg,'=',pos))
            g_tri_builder_replication_factor = parseInt (cfg,pos);

      	else if ((tok == "tri_accel_mb" || tok == "accel_mb") && parseSymbol (cfg,'=',pos))
            g_tri_accel = parseIdentifier (cfg,pos);
	else if ((tok == "tri_builder_mb" || tok == "builder_mb") && parseSymbol (cfg,'=',pos))
	    g_tri_builder = parseIdentifier (cfg,pos);
        else if ((tok == "tri_traverser_mb" || tok == "traverser_mb") && parseSymbol (cfg,'=',pos))
            g_tri_traverser = parseIdentifier (cfg,pos);

        else if (tok == "hair_accel" && parseSymbol (cfg,'=',pos))
            g_hair_accel = parseIdentifier (cfg,pos);
	else if (tok == "hair_builder" && parseSymbol (cfg,'=',pos))
            g_hair_builder = parseIdentifier (cfg,pos);
	else if (tok == "hair_traverser" && parseSymbol (cfg,'=',pos))
            g_hair_traverser = parseIdentifier (cfg,pos);
	else if (tok == "hair_builder_replication_factor" && parseSymbol (cfg,'=',pos))
            g_hair_builder_replication_factor = parseInt (cfg,pos);

        else if (tok == "subdiv_accel" && parseSymbol (cfg,'=',pos))
            g_subdiv_accel = parseIdentifier (cfg,pos);
	
        else if (tok == "verbose" && parseSymbol (cfg,'=',pos))
            g_verbose = parseInt (cfg,pos);
	else if (tok == "benchmark" && parseSymbol (cfg,'=',pos))
            g_benchmark = parseInt (cfg,pos);

        else if (tok == "flags") {
          g_scene_flags = 0;
          if (parseSymbol (cfg,'=',pos)) {
            do {
              std::string flag = parseIdentifier (cfg,pos);
              if      (flag == "static" ) g_scene_flags |= RTC_SCENE_STATIC;
              else if (flag == "dynamic") g_scene_flags |= RTC_SCENE_DYNAMIC;
              else if (flag == "compact") g_scene_flags |= RTC_SCENE_COMPACT;
              else if (flag == "coherent") g_scene_flags |= RTC_SCENE_COHERENT;
              else if (flag == "incoherent") g_scene_flags |= RTC_SCENE_INCOHERENT;
              else if (flag == "high_quality") g_scene_flags |= RTC_SCENE_HIGH_QUALITY;
              else if (flag == "robust") g_scene_flags |= RTC_SCENE_ROBUST;
            } while (parseSymbol (cfg,',',pos));
          }
        }
	else if (tok == "memory_preallocation_factor" && parseSymbol (cfg,'=',pos))
	  {
	    g_memory_preallocation_factor = parseFloat (cfg,pos);
	    DBG_PRINT( g_memory_preallocation_factor );
	  }

        else if (tok == "regression" && parseSymbol (cfg,'=',pos)) {
          g_regression_testing = parseInt (cfg,pos);
        }
        
      } while (findNext (cfg,',',pos));
    }

    if (g_verbose >= 1)
    {
      std::cout << "Embree Ray Tracing Kernels " << __EMBREE_VERSION__ << " (" << __DATE__ << ")" << std::endl;
      std::cout << "  Compiler : " << getCompilerName() << std::endl;
      std::cout << "  Platform : " << getPlatformName() << std::endl;
      std::cout << "  CPU      : " << stringOfCPUFeatures(getCPUFeatures()) << std::endl;
      std::cout << "  Features : ";
#if defined(RTCORE_RAY_MASK)
      std::cout << "raymasks ";
#endif
#if defined (RTCORE_BACKFACE_CULLING)
      std::cout << "backfaceculling ";
#endif
#if defined(RTCORE_INTERSECTION_FILTER)
      std::cout << "intersection_filter ";
#endif
#if defined(RTCORE_BUFFER_STRIDE)
      std::cout << "bufferstride ";
#endif
      std::cout << std::endl;

#if defined (__MIC__)
#if defined(RTCORE_BUFFER_STRIDE)
      std::cout << "  WARNING: enabled 'bufferstride' support will lower BVH build performance" << std::endl;
#endif
#endif
    }

    /* CPU has to support at least SSE2 */
#if !defined (__MIC__)
    if (!has_feature(SSE2)) {
      g_mutex.unlock();
      process_error(RTC_UNSUPPORTED_CPU,"CPU does not support SSE2");
      g_mutex.lock();
      return;
    }
#endif

    g_error = createTls();
    g_error_function = NULL;

    init_globals();

#if !defined(__MIC__)
    BVH4Register();
#else
    BVH4iRegister();
    BVH4MBRegister();
    BVH4HairRegister();

#endif 
#if defined(__TARGET_AVX__)
    if (has_feature(AVX)) {
      BVH8Register();
    }
#endif
    
    InstanceIntersectorsRegister();

    if (g_verbose >= 2) 
      printSettings();
    
    TaskScheduler::create(g_numThreads);

    /* execute regression tests */
    if (g_regression_testing) 
    {
      TaskScheduler::EventSync event;
      TaskScheduler::Task task(&event,task_regression_testing,NULL,TaskScheduler::getNumThreads(),NULL,NULL,"regression_testing");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
      event.sync();
    }

    CATCH_END;
  }
  
  RTCORE_API void rtcExit() 
  {
    Lock<MutexSys> lock(g_mutex);
    TRACE(rtcExit);
    CATCH_BEGIN;
    if (!g_initialized) {
      return;
    }
    TaskScheduler::destroy();
    {
      Lock<MutexSys> lock(g_errors_mutex);
      for (size_t i=0; i<g_errors.size(); i++)
        delete g_errors[i];
      destroyTls(g_error);
      g_errors.clear();
    }
    Alloc::global.clear();
    g_error_function = NULL;
    g_initialized = false;
    CATCH_END;
  }

  RTCError* getThreadError() 
  {
    RTCError* stored_error = (RTCError*) getTls(g_error);
    if (stored_error == NULL) {
      Lock<MutexSys> lock(g_errors_mutex);
      stored_error = new RTCError(RTC_NO_ERROR);
      g_errors.push_back(stored_error);
      setTls(g_error,stored_error);
    }
    return stored_error;
  }

  void process_error(RTCError error, const char* str)
  { 
    /* print error when in verbose mode */
    if (g_verbose) 
    {
      switch (error) {
      case RTC_NO_ERROR         : std::cerr << "Embree: No error"; break;
      case RTC_UNKNOWN_ERROR    : std::cerr << "Embree: Unknown error"; break;
      case RTC_INVALID_ARGUMENT : std::cerr << "Embree: Invalid argument"; break;
      case RTC_INVALID_OPERATION: std::cerr << "Embree: Invalid operation"; break;
      case RTC_OUT_OF_MEMORY    : std::cerr << "Embree: Out of memory"; break;
      case RTC_UNSUPPORTED_CPU  : std::cerr << "Embree: Unsupported CPU"; break;
      default                   : std::cerr << "Embree: Invalid error code"; break;                   
      };
      if (str) std::cerr << ", (" << str << ")";
      std::cerr << std::endl;
    }

    /* call user specified error callback */
    if (g_error_function) 
      g_error_function(error,str); 

    /* record error code */
    RTCError* stored_error = getThreadError();
    if (*stored_error == RTC_NO_ERROR)
      *stored_error = error;
  }

  RTCORE_API RTCError rtcGetError() 
  {
    TRACE(rtcGetError);
    RTCError* stored_error = getThreadError();
    RTCError error = *stored_error;
    *stored_error = RTC_NO_ERROR;
    return error;
  }

  RTCORE_API void rtcSetErrorFunction(RTC_ERROR_FUNCTION func) {
    g_error_function = func;
  }

  RTCORE_API void rtcDebug()
  {
    Lock<MutexSys> lock(g_mutex);

    TRACE(rtcDebug);
#if defined(RTCORE_STAT_COUNTERS)
    Stat::print(std::cout);
    Stat::clear();
#endif
#if defined(DEBUG) && 0
    extern void printTessCacheStats();
    printTessCacheStats();
#endif
  }
  
  RTCORE_API RTCScene rtcNewScene (RTCSceneFlags flags, RTCAlgorithmFlags aflags) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewScene);
    if (!isCoherent(flags) && !isIncoherent(flags)) flags = RTCSceneFlags(flags | RTC_SCENE_INCOHERENT);
    return (RTCScene) new Scene(flags,aflags);
    CATCH_END;
    return NULL;
  }
  
  RTCORE_API void rtcCommit (RTCScene scene) 
  {
    CATCH_BEGIN;
    TRACE(rtcCommit);
    VERIFY_HANDLE(scene);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.dumpGeometry(scene);
#endif

    ((Scene*)scene)->build(0,0);
    CATCH_END;
  }

  RTCORE_API void rtcCommitThread(RTCScene scene, unsigned int threadID, unsigned int numThreads) 
  {
    CATCH_BEGIN;
    TRACE(rtcCommitMT);
    VERIFY_HANDLE(scene);
    if (unlikely(numThreads == 0)) 
      process_error(RTC_INVALID_OPERATION,"invalid number of threads specified");

#if defined(__MIC__)
    if (unlikely(numThreads % 4 != 0 && numThreads != 1)) 
      FATAL("MIC requires numThreads % 4 == 0 in rtcCommitThread");
#endif
    
    ((Scene*)scene)->build(threadID,numThreads);

    CATCH_END;
  }
  
  RTCORE_API void rtcIntersect (RTCScene scene, RTCRay& ray) 
  {
    TRACE(rtcIntersect);
    STAT3(normal.travs,1,1,1);
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)&ray) & 0x0F) process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 16 bytes");   
#endif

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay old_ray = ray;
#endif

    ((Scene*)scene)->intersect(ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay1Intersect(scene,old_ray,ray);
#endif
  }
  
  RTCORE_API void rtcIntersect4 (const void* valid, RTCScene scene, RTCRay4& ray) 
  {
    TRACE(rtcIntersect4);
#if !defined(__TARGET_SIMD4__)
    process_error(RTC_INVALID_OPERATION,"rtcIntersect4 not supported");    
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x0F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 16 bytes");   
    if (((size_t)&ray ) & 0x0F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 16 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<4; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(normal.travs,1,cnt,4);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay4 old_ray = ray;
#endif

    ((Scene*)scene)->intersect4(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay4Intersect(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcIntersect8 (const void* valid, RTCScene scene, RTCRay8& ray) 
  {
    TRACE(rtcIntersect8);
#if !defined(__TARGET_SIMD8__)
    process_error(RTC_INVALID_OPERATION,"rtcIntersect8 not supported");                                    
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x1F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 32 bytes");   
    if (((size_t)&ray ) & 0x1F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 32 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<8; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(normal.travs,1,cnt,8);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay8 old_ray = ray;
#endif

    ((Scene*)scene)->intersect8(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay8Intersect(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcIntersect16 (const void* valid, RTCScene scene, RTCRay16& ray) 
  {
    TRACE(rtcIntersect16);
#if !defined(__TARGET_SIMD16__)
    process_error(RTC_INVALID_OPERATION,"rtcIntersect16 not supported");
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x3F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 64 bytes");   
    if (((size_t)&ray ) & 0x3F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 64 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<16; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(normal.travs,1,cnt,16);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay16 old_ray = ray;
#endif

    ((Scene*)scene)->intersect16(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay16Intersect(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcOccluded (RTCScene scene, RTCRay& ray) 
  {
    TRACE(rtcOccluded);
    STAT3(shadow.travs,1,1,1);
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)&ray) & 0x0F) process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 16 bytes");   
#endif

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay old_ray = ray;
#endif

    ((Scene*)scene)->occluded(ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay1Occluded(scene,old_ray,ray);
#endif

  }
  
  RTCORE_API void rtcOccluded4 (const void* valid, RTCScene scene, RTCRay4& ray) 
  {
    TRACE(rtcOccluded4);
#if !defined(__TARGET_SIMD4__)
    process_error(RTC_INVALID_OPERATION,"rtcOccluded4 not supported");
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x0F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 16 bytes");   
    if (((size_t)&ray ) & 0x0F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 16 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<4; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(shadow.travs,1,cnt,4);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay4 old_ray = ray;
#endif

    ((Scene*)scene)->occluded4(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay4Occluded(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcOccluded8 (const void* valid, RTCScene scene, RTCRay8& ray) 
  {
    TRACE(rtcOccluded8);
#if !defined(__TARGET_SIMD8__)
    process_error(RTC_INVALID_OPERATION,"rtcOccluded8 not supported");
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x1F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 32 bytes");   
    if (((size_t)&ray ) & 0x1F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 32 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<8; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(shadow.travs,1,cnt,8);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay8 old_ray = ray;
#endif

    ((Scene*)scene)->occluded8(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RayStreamLogger::rayStreamLogger.logRay8Occluded(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcOccluded16 (const void* valid, RTCScene scene, RTCRay16& ray) 
  {
    TRACE(rtcOccluded16);
#if !defined(__TARGET_SIMD16__)
    process_error(RTC_INVALID_OPERATION,"rtcOccluded16 not supported");
#else
#if defined(DEBUG)
    if (!((Scene*)scene)->is_build) process_error(RTC_INVALID_OPERATION,"scene got not committed");
    if (((size_t)valid) & 0x3F)  process_error(RTC_INVALID_ARGUMENT,"mask not aligned to 64 bytes");   
    if (((size_t)&ray ) & 0x3F)  process_error(RTC_INVALID_ARGUMENT,"ray not aligned to 64 bytes");   
#endif
    STAT(size_t cnt=0; for (size_t i=0; i<16; i++) cnt += ((int*)valid)[i] == -1;);
    STAT3(shadow.travs,1,cnt,16);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
    RTCRay16 old_ray = ray;
#endif

    ((Scene*)scene)->occluded16(valid,ray);

#if defined(RTCORE_ENABLE_RAYSTREAM_LOGGER)
  RayStreamLogger::rayStreamLogger.logRay16Occluded(valid,scene,old_ray,ray);
#endif

#endif
  }
  
  RTCORE_API void rtcDeleteScene (RTCScene scene) 
  {
    CATCH_BEGIN;
    TRACE(rtcDeleteScene);
    VERIFY_HANDLE(scene);
    delete (Scene*) scene;
    CATCH_END;
  }

  RTCORE_API unsigned rtcNewInstance (RTCScene target, RTCScene source) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewInstance);
    VERIFY_HANDLE(target);
    VERIFY_HANDLE(source);
    return ((Scene*) target)->newInstance((Scene*) source);
    CATCH_END;
    return -1;
  }

  RTCORE_API void rtcSetTransform (RTCScene scene, unsigned geomID, RTCMatrixType layout, const float* xfm) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetTransform);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    VERIFY_HANDLE(xfm);

    AffineSpace3fa transform = one;
    switch (layout) 
    {
    case RTC_MATRIX_ROW_MAJOR:
      transform = AffineSpace3fa(Vec3fa(xfm[ 0],xfm[ 4],xfm[ 8]),
                                 Vec3fa(xfm[ 1],xfm[ 5],xfm[ 9]),
                                 Vec3fa(xfm[ 2],xfm[ 6],xfm[10]),
                                 Vec3fa(xfm[ 3],xfm[ 7],xfm[11]));
      break;

    case RTC_MATRIX_COLUMN_MAJOR:
      transform = AffineSpace3fa(Vec3fa(xfm[ 0],xfm[ 1],xfm[ 2]),
                                 Vec3fa(xfm[ 3],xfm[ 4],xfm[ 5]),
                                 Vec3fa(xfm[ 6],xfm[ 7],xfm[ 8]),
                                 Vec3fa(xfm[ 9],xfm[10],xfm[11]));
      break;

    case RTC_MATRIX_COLUMN_MAJOR_ALIGNED16:
      transform = AffineSpace3fa(Vec3fa(xfm[ 0],xfm[ 1],xfm[ 2]),
                                 Vec3fa(xfm[ 4],xfm[ 5],xfm[ 6]),
                                 Vec3fa(xfm[ 8],xfm[ 9],xfm[10]),
                                 Vec3fa(xfm[12],xfm[13],xfm[14]));
      break;

    default: 
      process_error(RTC_INVALID_OPERATION,"Unknown matrix type");
      break;
    }
    ((Scene*) scene)->get_locked(geomID)->setTransform(transform);

    CATCH_END;
  }

  RTCORE_API unsigned rtcNewUserGeometry (RTCScene scene, size_t numItems) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewUserGeometry);
    VERIFY_HANDLE(scene);
    return ((Scene*)scene)->newUserGeometry(numItems);
    CATCH_END;
    return -1;
  }

  RTCORE_API unsigned rtcNewTriangleMesh (RTCScene scene, RTCGeometryFlags flags, size_t numTriangles, size_t numVertices, size_t numTimeSteps) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewTriangleMesh);
    VERIFY_HANDLE(scene);
    return ((Scene*)scene)->newTriangleMesh(flags,numTriangles,numVertices,numTimeSteps);
    CATCH_END;
    return -1;
  }

  RTCORE_API unsigned rtcNewHairGeometry (RTCScene scene, RTCGeometryFlags flags, size_t numCurves, size_t numVertices, size_t numTimeSteps) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewHairGeometry);
    VERIFY_HANDLE(scene);
    return ((Scene*)scene)->newBezierCurves(flags,numCurves,numVertices,numTimeSteps);
    CATCH_END;
    return -1;
  }

  RTCORE_API void rtcSetMask (RTCScene scene, unsigned geomID, int mask) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetMask);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setMask(mask);
    CATCH_END;
  }

  RTCORE_API void* rtcMapBuffer(RTCScene scene, unsigned geomID, RTCBufferType type) 
  {
    CATCH_BEGIN;
    TRACE(rtcMapBuffer);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    return ((Scene*)scene)->get_locked(geomID)->map(type);
    CATCH_END;
    return NULL;
  }

  RTCORE_API void rtcUnmapBuffer(RTCScene scene, unsigned geomID, RTCBufferType type) 
  {
    CATCH_BEGIN;
    TRACE(rtcUnmapBuffer);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->unmap(type);
    CATCH_END;
  }

  RTCORE_API void rtcSetBuffer(RTCScene scene, unsigned geomID, RTCBufferType type, void* ptr, size_t offset, size_t stride)
  {
    CATCH_BEGIN;
    TRACE(rtcSetBuffer);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setBuffer(type,ptr,offset,stride);
    CATCH_END;
  }

  RTCORE_API void rtcEnable (RTCScene scene, unsigned geomID) 
  {
    CATCH_BEGIN;
    TRACE(rtcEnable);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->enable();
    CATCH_END;
  }

  RTCORE_API void rtcUpdate (RTCScene scene, unsigned geomID) 
  {
    CATCH_BEGIN;
    TRACE(rtcUpdate);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->update();
    CATCH_END;
  }

  RTCORE_API void rtcUpdateBuffer (RTCScene scene, unsigned geomID, RTCBufferType type) 
  {
    CATCH_BEGIN;
    TRACE(rtcUpdateBuffer);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->updateBuffer(type);
    CATCH_END;
  }

  RTCORE_API void rtcDisable (RTCScene scene, unsigned geomID) 
  {
    CATCH_BEGIN;
    TRACE(rtcDisable);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->disable();
    CATCH_END;
  }

  RTCORE_API void rtcDeleteGeometry (RTCScene scene, unsigned geomID) 
  {
    CATCH_BEGIN;
    TRACE(rtcDeleteGeometry);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->erase();
    CATCH_END;
  }

  RTCORE_API void rtcSetUserData (RTCScene scene, unsigned geomID, void* ptr) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetUserData);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setUserData(ptr);
    CATCH_END;
  }

  RTCORE_API void rtcSetBoundsFunction (RTCScene scene, unsigned geomID, RTCBoundsFunc bounds)
  {
    CATCH_BEGIN;
    TRACE(rtcSetBoundsFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setBoundsFunction(bounds);
    CATCH_END;
  }

  RTCORE_API void rtcSetDisplacementFunction (RTCScene scene, unsigned geomID, RTCDisplacementFunc func, RTCBounds* bounds)
  {
    CATCH_BEGIN;
    TRACE(rtcSetDisplacementFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setDisplacementFunction(func,bounds);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectFunction (RTCScene scene, unsigned geomID, RTCIntersectFunc intersect) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectFunction(intersect);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectFunction4 (RTCScene scene, unsigned geomID, RTCIntersectFunc4 intersect4) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectFunction4);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectFunction4(intersect4);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectFunction8 (RTCScene scene, unsigned geomID, RTCIntersectFunc8 intersect8) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectFunction8);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectFunction8(intersect8);
    CATCH_END;
  }
  
  RTCORE_API void rtcSetIntersectFunction16 (RTCScene scene, unsigned geomID, RTCIntersectFunc16 intersect16) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectFunction16);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectFunction16(intersect16);
    CATCH_END;
  }

  RTCORE_API void rtcSetOccludedFunction (RTCScene scene, unsigned geomID, RTCOccludedFunc occluded) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOccludedFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOccludedFunction(occluded);
    CATCH_END;
  }

  RTCORE_API void rtcSetOccludedFunction4 (RTCScene scene, unsigned geomID, RTCOccludedFunc4 occluded4) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOccludedFunction4);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOccludedFunction4(occluded4);
    CATCH_END;
  }

  RTCORE_API void rtcSetOccludedFunction8 (RTCScene scene, unsigned geomID, RTCOccludedFunc8 occluded8) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOccludedFunction8);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOccludedFunction8(occluded8);
    CATCH_END;
  }

  RTCORE_API void rtcSetOccludedFunction16 (RTCScene scene, unsigned geomID, RTCOccludedFunc16 occluded16) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOccludedFunction16);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOccludedFunction16(occluded16);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectionFilterFunction (RTCScene scene, unsigned geomID, RTCFilterFunc intersect) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectionFilterFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectionFilterFunction(intersect);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectionFilterFunction4 (RTCScene scene, unsigned geomID, RTCFilterFunc4 filter4) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectionFilterFunction4);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectionFilterFunction4(filter4);
    CATCH_END;
  }

  RTCORE_API void rtcSetIntersectionFilterFunction8 (RTCScene scene, unsigned geomID, RTCFilterFunc8 filter8) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectionFilterFunction8);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectionFilterFunction8(filter8);
    CATCH_END;
  }
  
  RTCORE_API void rtcSetIntersectionFilterFunction16 (RTCScene scene, unsigned geomID, RTCFilterFunc16 filter16) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetIntersectionFilterFunction16);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setIntersectionFilterFunction16(filter16);
    CATCH_END;
  }

  RTCORE_API void rtcSetOcclusionFilterFunction (RTCScene scene, unsigned geomID, RTCFilterFunc intersect) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOcclusionFilterFunction);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOcclusionFilterFunction(intersect);
    CATCH_END;
  }

  RTCORE_API void rtcSetOcclusionFilterFunction4 (RTCScene scene, unsigned geomID, RTCFilterFunc4 filter4) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOcclusionFilterFunction4);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOcclusionFilterFunction4(filter4);
    CATCH_END;
  }

  RTCORE_API void rtcSetOcclusionFilterFunction8 (RTCScene scene, unsigned geomID, RTCFilterFunc8 filter8) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOcclusionFilterFunction8);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOcclusionFilterFunction8(filter8);
    CATCH_END;
  }
  
  RTCORE_API void rtcSetOcclusionFilterFunction16 (RTCScene scene, unsigned geomID, RTCFilterFunc16 filter16) 
  {
    CATCH_BEGIN;
    TRACE(rtcSetOcclusionFilterFunction16);
    VERIFY_HANDLE(scene);
    VERIFY_GEOMID(geomID);
    ((Scene*)scene)->get_locked(geomID)->setOcclusionFilterFunction16(filter16);
    CATCH_END;
  }

  /* new support for subdivision surfaces */
  RTCORE_API unsigned rtcNewSubdivisionMesh (RTCScene scene, RTCGeometryFlags flags, size_t numFaces, size_t numEdges, size_t numVertices, 
                                             size_t numEdgeCreases, size_t numVertexCreases, size_t numHoles, size_t numTimeSteps) 
  {
    CATCH_BEGIN;
    TRACE(rtcNewSubdivisionMesh);
    VERIFY_HANDLE(scene);
    return ((Scene*)scene)->newSubdivisionMesh(flags,numFaces,numEdges,numVertices,numEdgeCreases,numVertexCreases,numHoles,numTimeSteps);
    CATCH_END;
    return -1;
  }

}
