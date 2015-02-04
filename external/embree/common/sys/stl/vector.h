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

#include <stdio.h>
#include <assert.h>

#include "../platform.h"
#include "../ref.h" 

namespace embree
{
  template<class T>
    class vector_t : public RefCount
    {
    public:

      vector_t () : m_size(0), alloced(0), t(NULL) {}

      void clear() {
        if (t) alignedFree(t);
        m_size = alloced = 0;
        t = NULL;
      };

      vector_t(size_t sz) {
        m_size = 0; alloced = 0; t = NULL;
        if (sz) resize(sz);
      }

      vector_t(const vector_t<T> &other)
      {
        m_size = other.m_size;
        alloced = other.alloced;
        t = (T*)alignedMalloc(alloced*sizeof(T),64);
        for (size_t i=0; i<m_size; i++) t[i] = other.t[i];
      }
      
      ~vector_t() {
        if (t) alignedFree(t); t = NULL;
      }

      inline bool empty() const { return m_size == 0; }
      inline size_t size() const { return m_size; };

      T* begin() const { return t; };
      T* end() const { return t+m_size; };

	  __forceinline       T* data()       { return t; };
	  __forceinline const T* data() const { return t; };


      inline T& front() const { return t[0]; };
      inline T& back () const { return t[m_size-1]; };

      void push_back(const T &nt) {
        T v = nt; // need local copy as input reference could point to this vector
        reserve(m_size+1);
        t[m_size] = v;
        m_size++;
      }

	  void pop_back() {
        m_size--;
      }

      vector_t<T> &operator=(const vector_t<T> &other) {
        resize(other.m_size);
        for (size_t i=0;i<m_size;i++) t[i] = other.t[i];
        return *this;
      }

      __forceinline T& operator[](size_t i) {
        assert(t);
        assert(i < m_size);
        return t[i];
      };

      __forceinline const T& operator[](size_t i) const {
        assert(t);
        assert(i < m_size);
        return t[i];
      };

      void resize(size_t new_sz, bool exact = false)
      {
        if (new_sz < m_size) {
          if (exact) {
            T *old_t = t;
            t = (T*)alignedMalloc(new_sz*sizeof(T),64);
            for (size_t i=0;i<new_sz;i++) t[i] = old_t[i];
            alloced = new_sz;
            if (old_t) alignedFree(old_t);
          }
        } else {
          reserve(new_sz,exact);
        }
        m_size = new_sz;
      };

      void reserve(size_t sz, bool exact = false)
      {
        if (sz <= alloced) return;

        size_t newAlloced = alloced;
        if (exact) newAlloced = sz;
        else
          while (newAlloced < sz)
            newAlloced = (1 < (newAlloced * 2)) ? newAlloced * 2 : 1;

        T* old_t = t;
        assert(newAlloced > 0);
        t = (T*)alignedMalloc(newAlloced*sizeof(T),64);
        alloced = newAlloced;

        for (size_t i=0;i<m_size;i++) t[i] = old_t[i];

        if (old_t) alignedFree(old_t);
      }

    public:
      size_t m_size;    // number of valid items
      size_t alloced;   // number of items allocated
      T *t;             // data array
    };
}
