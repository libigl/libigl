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

#include "common/default.h"

namespace embree
{
  /*! Implements a data buffer. */
  class Buffer
  {
  public:

    /*! Buffer construction */
    Buffer (); 
    
    /*! Buffer destruction */
    ~Buffer ();
      
    /*! disallow copy */
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

  public:
    
    /*! initialized the buffer */
    void init(size_t num_in, size_t stride_in);

    /*! sets shared buffer */
    void set(void* ptr_in, size_t ofs_in, size_t stride_in);

    /*! sets shared buffer */
    void set(void* ptr);

    /*! allocated buffer */
    void alloc();

    /*! fees the buffer */
    void free();

    /*! maps the buffer */
    void* map(atomic_t& cntr);
    
    /*! unmaps the buffer */
    void unmap(atomic_t& cntr);

    /*! checks if the buffer is mapped */
    __forceinline bool isMapped() const {
      return mapped; 
    }

    /*! mark buffer as modified or unmodified */
    __forceinline void setModified(bool b) {
      modified = b;
    }

    /*! mark buffer as modified or unmodified */
    __forceinline bool isModified() const {
      return modified;
    }

    /*! returns the number of elements of the buffer */
    __forceinline size_t size() const { 
      return num; 
    }

    /*! returns true of the buffer is not empty */
    __forceinline operator bool() { 
      return ptr; 
    }

  protected:
    char* ptr;       //!< pointer to buffer data
    size_t bytes;    //!< size of buffer in bytes
    char* ptr_ofs;   //!< base pointer plus offset

#if defined(__MIC__)
    unsigned int stride;
#else
    size_t stride;   //!< stride of the stream in bytes
#endif
    size_t num;      //!< number of elements in the stream
    bool shared;     //!< set if memory is shared with application
    bool mapped;     //!< set if buffer is mapped
    bool modified;   //!< true if the buffer got modified
  };

  /*! Implements a data stream inside a data buffer. */
  template<typename T>
    class BufferT : public Buffer
  {
  public:

    typedef T value_type;

    /*! access to the ith element of the buffer stream */
    __forceinline const T& operator[](size_t i) const 
    {
      assert(i<num);
#if defined(RTCORE_BUFFER_STRIDE)
      return *(T*)(ptr_ofs + i*stride);
#else
      return *(T*)(ptr_ofs + i*sizeof(T));
#endif
    }

    __forceinline char* getPtr( size_t i = 0 ) const 
    {
#if defined(RTCORE_BUFFER_STRIDE)
      return ptr_ofs + i*stride;
#else
      return ptr_ofs + i*sizeof(T);
#endif
    }

    __forceinline unsigned int getBufferStride() const {
#if defined(RTCORE_BUFFER_STRIDE)
      return stride;
#else
      return sizeof(T);
#endif
    }
  };

  /*! Implements a data stream inside a data buffer. */
  template<>
    class BufferT<Vec3fa> : public Buffer
  {
  public:

    typedef Vec3fa value_type;

    /*! access to the ith element of the buffer stream */
    __forceinline const Vec3fa operator[](size_t i) const
    {
      assert(i<num);
#if defined(RTCORE_BUFFER_STRIDE)
#if defined(__MIC__)
      return *(Vec3fa*)(ptr_ofs + i*stride);
#else
      return Vec3fa(ssef::loadu(ptr_ofs + i*stride));
#endif
#else
      return *(Vec3fa*)(ptr_ofs + i*sizeof(Vec3fa));
#endif
    }

    __forceinline char* getPtr( size_t i = 0 ) const 
    {
#if defined(RTCORE_BUFFER_STRIDE)
      return ptr_ofs + i*stride;
#else
      return ptr_ofs + i*sizeof(Vec3fa);
#endif
    }

    __forceinline unsigned int getBufferStride() const {
#if defined(RTCORE_BUFFER_STRIDE)
      return stride;
#else
      return sizeof(Vec3fa);
#endif
    }
  };
}
