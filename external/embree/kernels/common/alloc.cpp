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

#include "alloc.h"

namespace embree
{
  Alloc Alloc::global;

  Alloc::Alloc () {
    blocks.reserve(1024);
  }

  Alloc::~Alloc () {
  }

  size_t Alloc::size() const {
    return size_t(blockSize)*blocks.size();
  }

  void Alloc::clear()
  {
    Lock<MutexSys> lock(mutex);
    for (size_t i=0; i<blocks.size(); i++)
      alignedFree(blocks[i]);
    blocks.clear();
  }
  
  void* Alloc::malloc() 
  {
    Lock<MutexSys> lock(mutex);
    if (blocks.size()) {
      void* ptr = blocks.back();
      blocks.pop_back();
      return ptr;
    }
    return alignedMalloc(blockSize,64);
  }
  
  void Alloc::free(void* ptr) 
  {
    Lock<MutexSys> lock(mutex);
    blocks.push_back(ptr);
  }
}

