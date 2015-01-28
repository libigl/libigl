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

#ifndef __EMBREE_PRIM_REF_H__
#define __EMBREE_PRIM_REF_H__

#include "common/default.h"

namespace embree
{
  /*! A primitive reference stores the bounds of the primitive and its ID. */
  struct __aligned(32) PrimRef 
  {
    __forceinline PrimRef () {}
    __forceinline PrimRef (const BBox3f& bounds, unsigned geomID, unsigned primID) {
      lower = bounds.lower; lower.a = geomID;
      upper = bounds.upper; upper.a = primID;
    }
    
    __forceinline const BBox3f bounds() const {
      return BBox3f(lower,upper);
    }

    __forceinline size_t geomID() const { 
      return lower.a;
    }

    __forceinline size_t primID() const { 
      return upper.a;
    }
    
#if defined(__MIC__)
    __forceinline void operator=(const PrimRef& v) { 
      const mic_f p = uload16f_low((float*)&v.lower);
      compactustore16f_low(0xff,(float*)this,p);            
    };
#endif

  public:
    Vec3fa lower;
    Vec3fa upper;
  };

  /*! Outputs primitive reference to a stream. */
  inline std::ostream& operator<<(std::ostream& cout, const PrimRef& ref) {
    return cout << "{ lower = " << ref.lower << ", upper = " << ref.upper << ", geomID = " << ref.geomID() << ", primID = " << ref.primID() << " }";
  }


  __forceinline void xchg(PrimRef& a, PrimRef& b)
  {
#if defined(__AVX__) || defined(__AVX2__)

    const avxf aa = load8f((float*)&a);
    const avxf bb = load8f((float*)&b);
    store8f((float*)&a,bb);
    store8f((float*)&b,aa);
#elif defined(__MIC__)
    const mic_f aa = uload16f_low((float*)&a.lower);
    const mic_f bb = uload16f_low((float*)&b.lower);
    compactustore16f_low(0xff,(float*)&b.lower,aa);
    compactustore16f_low(0xff,(float*)&a.lower,bb);
#else
    std::swap(a,b);
#endif
  }

  __forceinline bool subset(const PrimRef& a, const PrimRef& b)
  { 
    for ( size_t i = 0 ; i < 3 ; i++ ) if ( a.lower[i] < b.lower[i] ) return false;
    for ( size_t i = 0 ; i < 3 ; i++ ) if ( a.upper[i] > b.upper[i] ) return false;
    return true; 
  }


  __forceinline float area( const PrimRef& a ) 
  { 
    const Vec3fa d = a.upper - a.lower; 
    return 2.0f*(d.x*(d.y+d.z)+d.y*d.z); 
  }

}
#endif

