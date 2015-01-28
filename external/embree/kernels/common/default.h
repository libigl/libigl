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

#ifndef __EMBREE_RTCORE_XEON_DEFAULT_H__
#define __EMBREE_RTCORE_XEON_DEFAULT_H__

#include "sys/platform.h"
#include "sys/ref.h"
#include "sys/intrinsics.h"
#include "sys/sysinfo.h"
#include "sys/sync/atomic.h"
#include "sys/stl/vector.h"
#include "sys/stl/string.h"
#include "sys/taskscheduler.h"

#include "math/math.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include "math/bbox.h"
#include "math/affinespace.h"

#include "simd/simd.h"
#include "embree2/rtcore.h"
#include "stat.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>

#define __EMBREE_VERSION__ "2.2"

#define ERROR(x) \
  throw std::runtime_error(x)

#if defined(__EXIT_ON_ERROR__)
#define VERBOSE 1
#else
#define VERBOSE g_verbose
#endif

namespace embree
{
  /* global settings */
  extern size_t g_numThreads;
  extern size_t g_verbose;
  extern std::string g_top_accel;
  extern std::string g_tri_accel;
  extern std::string g_builder;
  extern std::string g_traverser;
  extern int g_scene_flags;
  extern size_t g_benchmark;

  /*! records an error */
  void recordError(RTCError error);

  /*! decoding of geometry flags */
  __forceinline bool isStatic    (RTCSceneFlags flags) { return (flags & 1) == RTC_SCENE_STATIC; }
  __forceinline bool isDynamic   (RTCSceneFlags flags) { return (flags & 1) == RTC_SCENE_DYNAMIC; }

  __forceinline bool isCompact   (RTCSceneFlags flags) { return flags & RTC_SCENE_COMPACT; }
  __forceinline bool isRobust    (RTCSceneFlags flags) { return flags & RTC_SCENE_ROBUST; }
  __forceinline bool isCoherent  (RTCSceneFlags flags) { return flags & RTC_SCENE_COHERENT; }
  __forceinline bool isIncoherent(RTCSceneFlags flags) { return flags & RTC_SCENE_INCOHERENT; }
  __forceinline bool isHighQuality(RTCSceneFlags flags) { return flags & RTC_SCENE_HIGH_QUALITY; }

  /*! CPU features */
  static const int SSE   = CPU_FEATURE_SSE; 
  static const int SSE2  = SSE | CPU_FEATURE_SSE2;
  static const int SSE3  = SSE2 | CPU_FEATURE_SSE3;
  static const int SSSE3 = SSE3 | CPU_FEATURE_SSSE3;
  static const int SSE41 = SSSE3 | CPU_FEATURE_SSE41;
  static const int SSE42 = SSE41 | CPU_FEATURE_SSE42 | CPU_FEATURE_POPCNT;
  static const int AVX   = SSE42 | CPU_FEATURE_AVX;
  static const int AVXI  = AVX | CPU_FEATURE_F16C | CPU_FEATURE_RDRAND;
  static const int AVX2  = AVXI | CPU_FEATURE_AVX2 | CPU_FEATURE_FMA3 | CPU_FEATURE_BMI1 | CPU_FEATURE_BMI2 | CPU_FEATURE_LZCNT;
  static const int KNC   = CPU_FEATURE_KNC;
  
  __forceinline bool has_feature(const int feature) {
    int cpu_features = getCPUFeatures();
    return (cpu_features & feature) == feature;
  }

#if defined (__MIC__)
#  define ISA KNC
#elif defined (__AVX2__)
#  define ISA AVX2
#elif defined(__AVXI__)
#  define ISA AVXI
#elif defined(__AVX__)
#  define ISA AVX
#elif defined (__SSE4_2__)
#  define ISA SSE42
#elif defined (__SSE4_1__)
#  define ISA SSE41
#elif defined(__SSSE3__)
#  define ISA SSSE3
#elif defined(__SSE3__)
#  define ISA SSE3
#elif defined(__SSE2__)
#  define ISA SSE2
#elif defined(__SSE__)
#  define ISA SSE
#elif defined (__MACOSX__)
#  define ISA SSSE3
#elif defined (__LINUX__)
#  define ISA SSE2
#elif defined (__WIN32__)
#  define ISA SSE2
#else 
#  define ISA SSE2
#endif

#if defined (__MACOSX__)
#if defined (__INTEL_COMPILER)
#define DEFAULT_ISA SSSE3
#else
#define DEFAULT_ISA SSE3
#endif
#else
#define DEFAULT_ISA SSE2
#endif

  inline std::string stringOfISA(int features)
  {
    if (features == SSE) return "SSE";
    if (features == SSE2) return "SSE2";
    if (features == SSE3) return "SSE3";
    if (features == SSSE3) return "SSSE3";
    if (features == SSE41) return "SSE4_1";
    if (features == SSE42) return "SSE4_2";
    if (features == AVX) return "AVX";
    if (features == AVXI) return "AVXI";
    if (features == AVX2) return "AVX2";
    if (features == KNC) return "KNC";
    return "UNKNOWN";
  }

#if defined (__SSE__) || defined (__MIC__)
  typedef Vec2<sseb> sse2b;
  typedef Vec3<sseb> sse3b;
  typedef Vec2<ssei> sse2i;
  typedef Vec3<ssei> sse3i;
  typedef Vec2<ssef> sse2f;
  typedef Vec3<ssef> sse3f;
#endif

#if defined (__AVX__)
  typedef Vec2<avxb> avx2b;
  typedef Vec3<avxb> avx3b;
  typedef Vec2<avxi> avx2i; 
  typedef Vec3<avxi> avx3i;
  typedef Vec2<avxf> avx2f;
  typedef Vec3<avxf> avx3f;
#endif

#if defined (__MIC__)
  typedef Vec2<mic_m> mic2b;
  typedef Vec3<mic_m> mic3b;
  typedef Vec2<mic_i> mic2i;
  typedef Vec3<mic_i> mic3i;
  typedef Vec2<mic_f> mic2f;
  typedef Vec3<mic_f> mic3f;
#endif

typedef void (*ErrorFunc) ();

#define DECLARE_SYMBOL(type,name)                  \
  namespace isa   { extern type name; }            \
  namespace sse41 { extern type name; }                                 \
  namespace avx   { extern type name; }                                 \
  namespace avx2  { extern type name; }                                 \
  void name##_error() { std::cerr << "Error: " << TOSTRING(name) << " not supported by your CPU" << std::endl; } \
  type name((type)name##_error);

#define SELECT_SYMBOL_DEFAULT(features,intersector) \
  intersector = isa::intersector;

#if defined(__TARGET_SSE41__)
#define SELECT_SYMBOL_SSE41(features,intersector) \
  if ((features & SSE41) == SSE41) intersector = sse41::intersector;
#else
#define SELECT_SYMBOL_SSE41(features,intersector)
#endif

#if defined(__TARGET_AVX__)
#define SELECT_SYMBOL_AVX(features,intersector) \
  if ((features & AVX) == AVX) intersector = avx::intersector;
#else
#define SELECT_SYMBOL_AVX(features,intersector)
#endif

#if defined(__TARGET_AVX2__)
#define SELECT_SYMBOL_AVX2(features,intersector) \
  if ((features & AVX2) == AVX2) intersector = avx2::intersector;
#else
#define SELECT_SYMBOL_AVX2(features,intersector)
#endif

#if defined(__MIC__)
#define SELECT_SYMBOL_KNC(features,intersector) \
  intersector = knc::intersector;
#else
#define SELECT_SYMBOL_KNC(features,intersector)
#endif

#define SELECT_SYMBOL_TEST(intersector)                               \
  if (!intersector) printf("WARNING: could not select code for " TOSTRING(intersector) "\n")
  
#define SELECT_SYMBOL_DEFAULT_SSE41(features,intersector) \
  SELECT_SYMBOL_DEFAULT(features,intersector);                                 \
  SELECT_SYMBOL_SSE41(features,intersector);                                   \
  SELECT_SYMBOL_TEST(intersector);

#define SELECT_SYMBOL_DEFAULT_AVX(features,intersector) \
  SELECT_SYMBOL_DEFAULT(features,intersector);                     \
  SELECT_SYMBOL_AVX(features,intersector);                         \
  SELECT_SYMBOL_TEST(intersector);

#define SELECT_SYMBOL_AVX_AVX2(features,intersector) \
  SELECT_SYMBOL_AVX(features,intersector);                         \
  SELECT_SYMBOL_AVX2(features,intersector);

#define SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,intersector) \
  SELECT_SYMBOL_DEFAULT(features,intersector);                     \
  SELECT_SYMBOL_AVX(features,intersector);                         \
  SELECT_SYMBOL_AVX2(features,intersector);                        \
  SELECT_SYMBOL_TEST(intersector);

#define SELECT_SYMBOL_DEFAULT_SSE41_AVX_AVX2(features,intersector) \
  SELECT_SYMBOL_DEFAULT(features,intersector);                     \
  SELECT_SYMBOL_SSE41(features,intersector);                       \
  SELECT_SYMBOL_AVX(features,intersector);                         \
  SELECT_SYMBOL_AVX2(features,intersector);                        \
  SELECT_SYMBOL_TEST(intersector);

#define SELECT_SYMBOL_DEFAULT_SSE41_AVX(features,intersector) \
  SELECT_SYMBOL_DEFAULT(features,intersector);                     \
  SELECT_SYMBOL_SSE41(features,intersector);                       \
  SELECT_SYMBOL_AVX(features,intersector);                         \
  SELECT_SYMBOL_TEST(intersector);
}

#endif
