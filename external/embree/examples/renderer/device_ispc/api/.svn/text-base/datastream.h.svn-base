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

#ifndef __EMBREE_ISPC_DATA_STREAM_H__
#define __EMBREE_ISPC_DATA_STREAM_H__

#include "../default.h"

namespace embree
{
  /*! A data stream. */
  class DataStream : public RefCount
  {
  public:

    DataStream (const Ref<Data>& ptr, size_t elements, size_t stride, size_t ofs) 
      : ptr(ptr), elements(elements), stride(stride), ofs(ofs) {}
    
    __forceinline size_t size() { return elements; }

    __forceinline Vec2f getVec2f(size_t i) {
      float* p = (float*)(ptr->map()+i*stride+ofs);
      return Vec2f(p[0],p[1]);
    }

    __forceinline Vector3f getVector3f(size_t i) {
      float* p = (float*)(ptr->map()+i*stride+ofs);
      return Vector3f(p[0],p[1],p[2]);
    }

    __forceinline Vector3i getVector3i(size_t i) {
      int* p = (int*)(ptr->map()+i*stride+ofs);
      return Vector3i(p[0],p[1],p[2]);
    }

  private:
    Ref<Data> ptr;    //!< Data buffer containing the stream.
    size_t elements;  //!< Number of elements in the stream.
    size_t stride;    //!< Stride in bytes between stream elements.
    size_t ofs;       //!< Start offset of first element in stream.
  };
}

#endif
