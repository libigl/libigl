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

#ifndef __EMBREE_AVX_H__
#define __EMBREE_AVX_H__

#include "simd/sse.h"

#if defined(__AVX__)
#include <immintrin.h>
#else
#include "immintrin_emu.h"
#endif

/* fix incompatibity issue of GCC and ICC immintrin.h */
#if defined(__INTEL_COMPILER) || defined(__EMU_M256_AVXIMMINTRIN_EMU_H__)
__forceinline __m256  _mm256_maskload_ps  (float const *ptr, __m256 mask) {
  return _mm256_maskload_ps(ptr, _mm256_castps_si256(mask));
}
__forceinline void    _mm256_maskstore_ps (float *ptr, __m256 mask, __m256 data) {
  _mm256_maskstore_ps(ptr, _mm256_castps_si256(mask), data);
}
#elif !defined(_MSC_VER) && !defined(__UNIX__)
__forceinline __m256  _mm256_maskload_ps  (float const *ptr, __m256i mask) {
  return _mm256_maskload_ps(ptr, _mm256_castsi256_ps(mask));
}
__forceinline void    _mm256_maskstore_ps (float *ptr, __m256i mask, __m256 data) {
  _mm256_maskstore_ps(ptr, _mm256_castsi256_ps(mask), data);
}
#endif

#include "simd/avxb.h"
#include "simd/avxi.h"
#include "simd/avxf.h"

#define BEGIN_ITERATE_AVXB(valid_i,id_o) { \
  int _valid_t = movemask(valid_i);                       \
  while (_valid_t) {                                      \
    int id_o = __bsf(_valid_t);                            \
    _valid_t = __btc(_valid_t,id_o);
#define END_ITERATE_AVXB } }

#define BEGIN_ITERATE_AVXI(valid_i,obj_i,valid_o,obj_o) { \
  int _valid_t = movemask(valid_i);                       \
  while (_valid_t) {                                      \
    int obj_o = obj_i[__bsf(_valid_t)];                         \
    avxb valid_o = valid_i & (obj_i == broadcast(&obj_i[__bsf(_valid_t)])); \
    _valid_t ^= movemask(valid_o);
#define END_ITERATE_AVXI } }

#endif
