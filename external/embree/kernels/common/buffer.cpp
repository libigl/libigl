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

#include "buffer.h"

namespace embree
{
  Buffer::Buffer () 
    : ptr(NULL), bytes(0), ptr_ofs(NULL), stride(0), num(0), shared(false), mapped(false), modified(true) {}
  
  Buffer::~Buffer () {
    free();
  }

  void Buffer::init(size_t num_in, size_t stride_in) 
  {
    ptr = NULL;
    bytes = num_in*stride_in;
    ptr_ofs = NULL;
    num = num_in;
    stride = stride_in;
    shared = false;
    mapped = false;
    modified = true;
  }

  void Buffer::set(void* ptr_in, size_t ofs_in, size_t stride_in)
  {
#if !defined(RTCORE_BUFFER_STRIDE)
    if (stride_in != stride) {
      process_error(RTC_INVALID_OPERATION,"buffer stride feature disabled at compile time and specified stride does not match default stride");
      return;
    }
#endif

    ptr = (char*) ptr_in;
    bytes = 0;
    ptr_ofs = (char*) ptr_in + ofs_in;
    stride = stride_in;
    shared = true;
  }

  void Buffer::alloc() {
    ptr = ptr_ofs = (char*) alignedMalloc(bytes);
  }

  void Buffer::free()
  {
    if (shared || !ptr) return;
    alignedFree(ptr); ptr = NULL; ptr_ofs = NULL; bytes = 0;
  }

  void* Buffer::map(atomic_t& cntr)
  {
    /* report error if buffer is already mapped */
    if (mapped) {
      process_error(RTC_INVALID_OPERATION,"buffer is already mapped");
      return NULL;
    }

    /* allocate buffer */
    if (!ptr && !shared && bytes)
      alloc();

    /* return mapped buffer */
    atomic_add(&cntr,+1); 
    mapped = true;
    return ptr;
  }

  void Buffer::unmap(atomic_t& cntr)
  {
    /* report error if buffer not mapped */
    if (!mapped) {
      process_error(RTC_INVALID_OPERATION,"buffer is not mapped");
      return;
    }

    /* unmap buffer */
    atomic_add(&cntr,-1); 
    mapped = false;
  }
}
