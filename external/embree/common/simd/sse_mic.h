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

#ifndef __EMBREE_SSE_MIC_H__
#define __EMBREE_SSE_MIC_H__

#include "sys/platform.h"
#include "sys/intrinsics.h"

#include <immintrin.h>

namespace embree {
  struct sseb_t;
  struct ssei_t;
  struct ssef_t;
  struct sseb_m;
  struct ssei_m;
  struct ssef_m;
}

#include "simd/sseb_mic.h"
#include "simd/ssei_mic.h"
#include "simd/ssef_mic.h"

namespace embree {
  typedef sseb_m sseb;
  typedef ssei_m ssei;
  typedef ssef_m ssef;

__forceinline const ssei cast( const ssef& a ) { return *(ssei*)&a; }
__forceinline const ssef cast( const ssei& a ) { return *(ssef*)&a; }

 __forceinline ssei_t floor_i( const ssef_t& other ) 
 { 
   const __m512i m512 = _mm512_cvtfxpnt_round_adjustps_epi32(other.m512,_MM_FROUND_FLOOR,_MM_EXPADJ_NONE); 
   return ssei_t(m512);
 }

}




#endif
