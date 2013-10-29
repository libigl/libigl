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

#include "../common/default.h"

namespace embree
{
  /*! A primitive reference stores the bounds of the primitive and its ID. */
  struct __align(16) PrimRef 
  {
    __forceinline PrimRef () {}
    __forceinline PrimRef (const BBox3f& bounds, size_t id) {
      assert(!bounds.empty() || id == size_t(-1));
      lower = bounds.lower; 
      upper = bounds.upper;
#if defined(__X86_64__)
      lower.a = (unsigned) (id >>  0);
      upper.a = (unsigned) (id >> 32);
#else
      lower.a = id;
#endif
    }
    
    __forceinline size_t id() const { 
#if defined(__X86_64__)
      return (size_t(unsigned(upper.a)) << 32) | size_t(unsigned(lower.a));
#else
      return lower.a;
#endif
    }
    
    __forceinline const BBox3f bounds() const {
      return BBox3f(lower,upper);
    }

  public:
    Vec3fa lower;
    Vec3fa upper;
  };

  /*! Outputs primitive reference to a stream. */
  inline std::ostream& operator<<(std::ostream& cout, const PrimRef& ref) {
    return cout << "{ lower = " << ref.lower << ", upper = " << ref.upper << ", id = " << ref.id() << " }";
  }
}
#endif

