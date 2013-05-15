// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_SSE_H__
#define __EMBREE_SSE_H__

#include "sys/platform.h"
#include "sys/intrinsics.h"

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSE3__
#include <pmmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

const __m128 _mm_lookupmask_ps[16] = {
  _mm_castsi128_ps(_mm_set_epi32( 0, 0, 0, 0)),
  _mm_castsi128_ps(_mm_set_epi32( 0, 0, 0,-1)),
  _mm_castsi128_ps(_mm_set_epi32( 0, 0,-1, 0)),
  _mm_castsi128_ps(_mm_set_epi32( 0, 0,-1,-1)),
  _mm_castsi128_ps(_mm_set_epi32( 0,-1, 0, 0)),
  _mm_castsi128_ps(_mm_set_epi32( 0,-1, 0,-1)),
  _mm_castsi128_ps(_mm_set_epi32( 0,-1,-1, 0)),
  _mm_castsi128_ps(_mm_set_epi32( 0,-1,-1,-1)),
  _mm_castsi128_ps(_mm_set_epi32(-1, 0, 0, 0)),
  _mm_castsi128_ps(_mm_set_epi32(-1, 0, 0,-1)),
  _mm_castsi128_ps(_mm_set_epi32(-1, 0,-1, 0)),
  _mm_castsi128_ps(_mm_set_epi32(-1, 0,-1,-1)),
  _mm_castsi128_ps(_mm_set_epi32(-1,-1, 0, 0)),
  _mm_castsi128_ps(_mm_set_epi32(-1,-1, 0,-1)),
  _mm_castsi128_ps(_mm_set_epi32(-1,-1,-1, 0)),
  _mm_castsi128_ps(_mm_set_epi32(-1,-1,-1,-1))
};

#if defined (__SSE4_1__) || defined (__SSE4_2__)
#include <smmintrin.h>
#else
#include "smmintrin_emu.h"
#endif

#if defined (__AES__) || defined (__PCLMUL__)
#include <wmmintrin.h>
#endif

#include "simd/sseb.h"
#include "simd/ssei.h"
#include "simd/ssef.h"

#define BEGIN_ITERATE_SSEB(valid_i,id_o) { \
  int _valid_t = movemask(valid_i);                       \
  while (_valid_t) {                                      \
    int id_o = __bsf(_valid_t);                            \
    _valid_t = __btc(_valid_t,id_o);
#define END_ITERATE_SSEB } }

#define BEGIN_ITERATE_SSEI(valid_i,obj_i,valid_o,obj_o) { \
  int _valid_t = movemask(valid_i);                       \
  while (_valid_t) {                                      \
    int obj_o = obj_i[__bsf(_valid_t)];                         \
    sseb valid_o = valid_i & (obj_i == ssei(obj_o));            \
    _valid_t ^= movemask(valid_o);
#define END_ITERATE_SSEI } }

namespace embree
{
  /* prefetches */
  __forceinline void prefetchL1 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T0); }
  __forceinline void prefetchL2 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T1); }
  __forceinline void prefetchL3 (const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_T2); }
  __forceinline void prefetchNTA(const void* ptr) { _mm_prefetch((const char*)ptr,_MM_HINT_NTA); }
}

#endif
