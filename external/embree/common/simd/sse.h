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

#pragma once

#include "sys/platform.h"
#include "sys/intrinsics.h"
#include "sse_special.h"

namespace embree 
{
#if defined(__SSE4_1__)
  __forceinline __m128 blendv_ps( __m128 f, __m128 t, __m128 mask ) { 
    return _mm_blendv_ps(f,t,mask);
  }
#else
  __forceinline __m128 blendv_ps( __m128 f, __m128 t, __m128 mask ) { 
    return _mm_or_ps(_mm_and_ps(mask, t), _mm_andnot_ps(mask, f)); 
  }
#endif

  extern const __m128 _mm_lookupmask_ps[16];

  struct sseb;
  struct ssei;
  struct ssef;

#if !defined(__MIC__)
  typedef ssef ssef_t;
  typedef ssei ssei_t;

  typedef ssef ssef_m;
  typedef ssei ssei_m;
#endif
}

#include "simd/sseb.h"
#include "simd/ssei.h"
#include "simd/ssef.h"
