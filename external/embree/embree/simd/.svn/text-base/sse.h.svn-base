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

#if defined (__SSE4_1__)
#include <smmintrin.h>

/*! workaround for compiler bug in VS2008 */
#if defined(_MSC_VER) && (_MSC_VER < 1600)
  #define _mm_blendv_ps __emu_mm_blendv_ps
  __forceinline __m128 _mm_blendv_ps( __m128 f, __m128 t, __m128 mask ) { 
    return _mm_or_ps(_mm_and_ps(mask, t), _mm_andnot_ps(mask, f)); 
  }
  #define _mm_extract_epi32 __emu_mm_extract_epi32
  __forceinline int _mm_extract_epi32( __m128i input, const int i ) {
    //return _mm_cvtsi128_si32(_mm_shuffle_epi32(input, _MM_SHUFFLE(i, i, i, i)));
    return input.m128i_i32[i];
  }
#endif

#else
#include "smmintrin_emu.h"
#endif

#if defined (__SSE4_2__)
#include <nmmintrin.h>
#endif

#if defined (__AES__)
#include <wmmintrin.h>
#endif

namespace embree 
{
  struct sseb;
  struct ssei;
  struct ssef;
}

#include "simd/sseb.h"
#include "simd/ssei.h"
#include "simd/ssef.h"

namespace embree 
{
  typedef sseb sseb_t;
  typedef ssei ssei_t;
  typedef ssef ssef_t;

  typedef sseb sseb_m;
  typedef ssei ssei_m;
  typedef ssef ssef_m;
}

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

#endif
