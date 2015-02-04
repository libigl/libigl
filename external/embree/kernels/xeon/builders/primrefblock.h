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

#include "common/primref.h"
#include "common/atomic_set.h"

namespace embree
{
  /*! Block of build primitives */
  template<typename PrimRef>
    class PrimRefBlockT
  {
  public:
    
    typedef PrimRef T;
    
    /*! Number of primitive references inside a block */
    static const size_t blockSize = 511;
    
    /*! default constructor */
    PrimRefBlockT () : num(0) {}
    
    /*! frees the block */
    __forceinline void clear(size_t n = 0) { num = n; }
    
    /*! return base pointer */
    __forceinline PrimRef* base() { return ptr; }
    
    /*! returns number of elements */
    __forceinline size_t size() const { return num; }
    
    /*! inserts a primitive reference */
    __forceinline bool insert(const PrimRef& ref) {
      if (unlikely(num >= blockSize)) return false;
      ptr[num++] = ref;
      return true;
    }
    
    /*! access the i-th primitive reference */
    __forceinline       PrimRef& operator[] (size_t i)       { return ptr[i]; }
    __forceinline const PrimRef& operator[] (size_t i) const { return ptr[i]; }
    
    /*! access the i-th primitive reference */
    __forceinline       PrimRef& at (size_t i)       { return ptr[i]; }
    __forceinline const PrimRef& at (size_t i) const { return ptr[i]; }
    
  private:
    PrimRef ptr[blockSize];   //!< Block with primitive references
    size_t num;               //!< Number of primitive references in block
  };
  
  typedef PrimRefBlockT<PrimRef> PrimRefBlock;
}
