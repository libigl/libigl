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

#ifndef __EMBREE_AVX_H__
#define __EMBREE_AVX_H__

#include "simd/sse.h"

#if defined(__AVX__)
#include <immintrin.h>
#else
#include "immintrin_emu.h"
#endif

#if defined (__AVX_I__)
#endif

#if defined (__AVX2__)
#endif

namespace embree 
{
  struct avxb;
  struct avxi;
  struct avxf;
}

#include "simd/avxb.h"
#if defined (__AVX_I__)
#include "simd/avxi.h"
#else
#include "simd/avxi_emu.h"
#endif
#include "simd/avxf.h"

namespace embree 
{
  typedef avxb avxb_t;
  typedef avxi avxi_t;
  typedef avxf avxf_t;

  typedef avxb avxb_m;
  typedef avxi avxi_m;
  typedef avxf avxf_m;
}

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
