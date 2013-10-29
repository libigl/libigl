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

#ifndef MIC_M_H
#define MIC_M_H

#include "sys/platform.h"

namespace embree
{
  class __align(4) mic_m
  {
  public:
    __mmask v;
    
    /* construction */
    __forceinline mic_m() {};
    __forceinline mic_m(bool b) { v = b ? 0xFFFF : 0x0000; };
    __forceinline mic_m(const __mmask &t) { v = t; };
    __forceinline mic_m(const mic_m &t) { v = t; };
    __forceinline mic_m(int t) { v = (__mmask)t; };
    __forceinline mic_m(unsigned int t) { v = (__mmask)t; };
    __forceinline mic_m(unsigned long t) { v = (__mmask)t; };
    
    __forceinline operator __mmask () const { return v; };
    //__forceinline operator unsigned () const { return _mm512_mask2int(v); };
    
    /* assignment operators */
    __forceinline void operator=(const mic_m &f) { v = f.v; };
    
    __forceinline void operator&=(const mic_m &a) { v = _mm512_kand(v,a); }; 
    __forceinline void operator|=(const mic_m &a) { v = _mm512_kor(v,a); }; 
    __forceinline void operator^=(const mic_m &a) { v = _mm512_kxor(v,a); }; 
  };  

  /* constants */
#define MIC_M_ZERO mic_m(0x0000)
#define MIC_M_ALL  mic_m(0xffff)
  
# define MIC_NO_BIT_SET_32 32
# define MIC_NO_BIT_SET_64 64
  
  
  __forceinline unsigned int toInt (const mic_m &a) { return _mm512_mask2int(a); };
  __forceinline mic_m        toMask(const int &a)   { return _mm512_int2mask(a); };
  
  /* comparisons */
  __forceinline mic_m bitwise_eq(const mic_m &a, const mic_m &b) { return _mm512_kxnor((__mmask)a,(__mmask)b); }
  __forceinline mic_m bitwise_ne(const mic_m &a, const mic_m &b) { return _mm512_kxor((__mmask)a,(__mmask)b); }
  
  __forceinline int nz_mask(const mic_m &a) { return !_mm512_kortestz(a,a); }
  __forceinline int ez_mask(const mic_m &a) { return  _mm512_kortestz(a,a); }
  
  __forceinline int eq_mask(const mic_m &a,const mic_m &b)
  {
    const mic_m c = _mm512_kxor((__mmask)a,(__mmask)b);
    return ez_mask(c);
  }
  
  __forceinline int neq_mask(const mic_m &a,const mic_m &b)
  {
    const mic_m c = _mm512_kxor((__mmask)a,(__mmask)b);
    return nz_mask(c);  
  }
  
  __forceinline mic_m andn(const mic_m &a, const mic_m &b) { return _mm512_kandn(b,a); }
  
  /* bit count */
  __forceinline unsigned long countbits(unsigned long r1) { 
    return _mm_countbits_64(r1); 
  };
  
  __forceinline unsigned int countbits(unsigned int r1) {
    return _mm_countbits_32(r1); 
  };
  
  __forceinline unsigned int countbits(mic_m r1) { 
    return _mm_countbits_32((unsigned int)r1); 
  };
  
  /* reductions */
  __forceinline bool all(const mic_m &a)  { return  _mm512_kortestc(a,a); }
  __forceinline bool any(const mic_m &a)  { return !_mm512_kortestz(a,a); }
  __forceinline bool none(const mic_m &a) { return  _mm512_kortestz(a,a); }
  
  /* selection */
  __forceinline mic_m sel(const mic_m &s, const mic_m &a, const mic_m &b) {
    return _mm512_kor(_mm512_kand((__mmask)s, (__mmask)a), _mm512_kandn((__mmask)s, (__mmask)b));
  }
  
  __forceinline unsigned long bsf64(const unsigned long mask) 
  { 
    return _mm_tzcnt_64(mask); 
  };
  
  __forceinline unsigned int bsf32(const unsigned int mask) 
  { 
    return _mm_tzcnt_32(mask); 
  };
  
  __forceinline unsigned long bsf64(const long start,const unsigned long mask) { 
    return _mm_tzcnti_64(start,mask); 
  };
  
  __forceinline mic_m operator&(const mic_m &a, const mic_m &b) { return _mm512_kand((__mmask)a,(__mmask)b); };
  __forceinline mic_m operator|(const mic_m &a, const mic_m &b) { return _mm512_kor((__mmask)a,(__mmask)b); };
  __forceinline mic_m operator^(const mic_m &a, const mic_m &b) { return _mm512_kxor((__mmask)a,(__mmask)b); };
  __forceinline mic_m operator~(const mic_m &a) { return _mm512_knot((__mmask)a); }
  
  __forceinline std::ostream& operator<<(std::ostream& o, const mic_m& v) {
    return o << "(" << (unsigned short)v << ")";
  }
}

#endif
