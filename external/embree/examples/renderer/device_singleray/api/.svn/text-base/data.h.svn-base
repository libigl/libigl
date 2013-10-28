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

#ifndef __EMBREE_DATA_H__
#define __EMBREE_DATA_H__

#include "../default.h"

namespace embree
{
  /*! Data container. */
  class Data : public RefCount
  {
  public:

    Data (size_t bytes) : bytes(bytes) {
      ptr = alignedMalloc(bytes);
    }
    
    Data (size_t bytes, const void* ptr_i, bool copy = true) : bytes(bytes) {
      if (copy) {
        ptr = alignedMalloc(bytes);
        memcpy(ptr,ptr_i,bytes);
      } else {
        ptr = (void*)ptr_i;
      }
    }
    
    virtual ~Data () {
      alignedFree(ptr); ptr = NULL;
      bytes = 0;
    }

    __forceinline const char* map() const { return (const char*)ptr; }
    __forceinline       char* map()       { return (      char*)ptr; }

    size_t size() { return bytes; }
    
  private:
    void* ptr;
    size_t bytes;
  };
}

#endif
