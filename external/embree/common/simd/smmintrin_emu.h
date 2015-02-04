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

#pragma once // FIXME: remove this file?

#ifndef __GNUC__
#pragma message (" --- Intel remark: SSE4 intrinsics are emulated with SSE3 ---")
#endif

#define _MM_FROUND_TO_NEAREST_INT    0x00
#define _MM_FROUND_TO_NEG_INF        0x01
#define _MM_FROUND_TO_POS_INF        0x02
#define _MM_FROUND_TO_ZERO           0x03
#define _MM_FROUND_CUR_DIRECTION     0x04

#define _mm_blendv_ps __emu_mm_blendv_ps
__forceinline __m128 _mm_blendv_ps( __m128 value, __m128 input, __m128 mask ) { 
    return _mm_or_ps(_mm_and_ps(mask, input), _mm_andnot_ps(mask, value)); 
}

#define _mm_blend_ps __emu_mm_blend_ps
__forceinline __m128 _mm_blend_ps( __m128 value, __m128 input, const int mask ) { 
    assert(mask < 0x10); return _mm_blendv_ps(value, input, _mm_lookupmask_ps[mask]); 
}

#define _mm_blendv_epi8 __emu_mm_blendv_epi8
__forceinline __m128i _mm_blendv_epi8( __m128i value, __m128i input, __m128i mask ) { 
    return _mm_or_si128(_mm_and_si128(mask, input), _mm_andnot_si128(mask, value)); 
}

#define _mm_mullo_epi32 __emu_mm_mullo_epi32
__forceinline __m128i _mm_mullo_epi32( __m128i value, __m128i input ) {
  __m128i rvalue;
  char* _r = (char*)(&rvalue + 1);
  char* _v = (char*)(& value + 1);
  char* _i = (char*)(& input + 1);
  for ( ssize_t i = -16 ; i != 0 ; i += 4 ) *((int32*)(_r + i)) = *((int32*)(_v + i))*  *((int32*)(_i + i));
  return rvalue;
}


#define _mm_min_epi32 __emu_mm_min_epi32
__forceinline __m128i _mm_min_epi32( __m128i value, __m128i input ) { 
    return _mm_blendv_epi8(input, value, _mm_cmplt_epi32(value, input)); 
}

#define _mm_max_epi32 __emu_mm_max_epi32
__forceinline __m128i _mm_max_epi32( __m128i value, __m128i input ) { 
    return _mm_blendv_epi8(value, input, _mm_cmplt_epi32(value, input)); 
}

#define _mm_extract_epi32 __emu_mm_extract_epi32
__forceinline int _mm_extract_epi32( __m128i input, const int index ) {
  switch ( index ) {
  case 0: return _mm_cvtsi128_si32(input);
  case 1: return _mm_cvtsi128_si32(_mm_shuffle_epi32(input, _MM_SHUFFLE(1, 1, 1, 1)));
  case 2: return _mm_cvtsi128_si32(_mm_shuffle_epi32(input, _MM_SHUFFLE(2, 2, 2, 2)));
  case 3: return _mm_cvtsi128_si32(_mm_shuffle_epi32(input, _MM_SHUFFLE(3, 3, 3, 3)));
  default: assert(false); return 0;
  }
}

#define _mm_insert_epi32 __emu_mm_insert_epi32
__forceinline __m128i _mm_insert_epi32( __m128i value, int input, const int index ) { 
    assert(index >= 0 && index < 4); ((int*)&value)[index] = input; return value; 
}

#define _mm_extract_ps __emu_mm_extract_ps
__forceinline int _mm_extract_ps( __m128 input, const int index ) {
  int32* ptr = (int32*)&input; return ptr[index];
}

#define _mm_insert_ps __emu_mm_insert_ps
__forceinline __m128 _mm_insert_ps( __m128 value, __m128 input, const int index )
{ assert(index < 0x100); ((float*)&value)[(index >> 4)&0x3] = ((float*)&input)[index >> 6]; return _mm_andnot_ps(_mm_lookupmask_ps[index&0xf], value); }

#define _mm_round_ps __emu_mm_round_ps
__forceinline __m128 _mm_round_ps( __m128 value, const int flags )
{
  switch ( flags )
  {
  case _MM_FROUND_TO_NEAREST_INT: return _mm_cvtepi32_ps(_mm_cvtps_epi32(value));
  case _MM_FROUND_TO_NEG_INF    : return _mm_cvtepi32_ps(_mm_cvtps_epi32(_mm_add_ps(value, _mm_set1_ps(-0.5f))));
  case _MM_FROUND_TO_POS_INF    : return _mm_cvtepi32_ps(_mm_cvtps_epi32(_mm_add_ps(value, _mm_set1_ps( 0.5f))));
  case _MM_FROUND_TO_ZERO       : return _mm_cvtepi32_ps(_mm_cvttps_epi32(value));
  }
  return value;
}

#ifdef _M_X64
#define _mm_insert_epi64 __emu_mm_insert_epi64
__forceinline __m128i _mm_insert_epi64( __m128i value, __int64 input, const int index ) { 
    assert(size_t(index) < 4); ((__int64*)&value)[index] = input; return value; 
}

#define _mm_extract_epi64 __emu_mm_extract_epi64
__forceinline __int64 _mm_extract_epi64( __m128i input, const int index ) { 
    assert(size_t(index) < 2); 
    return index == 0 ? _mm_cvtsi128_si64x(input) : _mm_cvtsi128_si64x(_mm_unpackhi_epi64(input, input)); 
}
#endif
