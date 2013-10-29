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

#ifndef __EMBREE_RTCORE_DEFAULT_H__
#define __EMBREE_RTCORE_DEFAULT_H__

#include "sys/platform.h"
#include "sys/ref.h"
#include "sys/intrinsics.h"
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

#include "stat.h"

#include <map>
#include <set>

namespace embree
{
#if defined (__SSE__)
  typedef Vec2<sseb> sse2b;
  typedef Vec3<sseb> sse3b;
  typedef Vec2<ssei> sse2i;
  typedef Vec3<ssei> sse3i;
  typedef Vec2<ssef> sse2f;
  typedef Vec3<ssef> sse3f;

  typedef Vec3<sseb_m> sse3b_m;
  typedef Vec3<ssei_m> sse3i_m;
  typedef Vec3<ssef_m> sse3f_m;
#endif

#if defined (__AVX__)
  typedef Vec2<avxb_t> avx2b;
  typedef Vec3<avxb_t> avx3b;
  typedef Vec2<avxi_t> avx2i; 
  typedef Vec3<avxi_t> avx3i;
  typedef Vec2<avxf_t> avx2f;
  typedef Vec3<avxf_t> avx3f;

  typedef Vec2<avxb_m> avx2b_m;
  typedef Vec3<avxb_m> avx3b_m;
  typedef Vec2<avxi_m> avx2i_m;
  typedef Vec3<avxi_m> avx3i_m;
  typedef Vec2<avxf_m> avx2f_m;
  typedef Vec3<avxf_m> avx3f_m;
#endif

#if defined(__MIC__)
  typedef Vec2<mic_m> mic2b;
  typedef Vec3<mic_m> mic3b;
  typedef Vec2<mic_i> mic2i;
  typedef Vec3<mic_i> mic3i;
  typedef Vec2<mic_f> mic2f;
  typedef Vec3<mic_f> mic3f;
#endif

  extern int g_verbose;
}

#endif
