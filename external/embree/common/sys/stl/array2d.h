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

#include "../platform.h"

namespace embree
{
  template<typename T>
    class Array2D
  {
  public:
    Array2D () : array(NULL), size_x(0), size_y(0) {}
    ~Array2D () { alignedFree(array); }
    
    void init(size_t width, size_t height) {
      size_x = width; size_y = height;
      //delete[] array; array = new T[width*height];
      alignedFree(array); array = (T*)alignedMalloc(width*height*sizeof(T));
    }
    
    void init(size_t width, size_t height, const T& v) 
    {
      size_x = width; size_y = height;
      //delete[] array; array = new T[width*height];
      alignedFree(array); array = (T*)alignedMalloc(width*height*sizeof(T));
      for (size_t y=0; y<height; y++)
        for (size_t x=0; x<width; x++)
          array[y*size_x+x] = v;
    }
    
    __forceinline size_t width() const { 
      return size_x;
    }
    
    __forceinline size_t height() const {
      return size_y;
    }
    
    __forceinline T& operator() (size_t x, size_t y) {
      assert(x<size_x);
      assert(y<size_y);
      return array[y*size_x+x];
    }

    __forceinline const T& operator() (size_t x, size_t y) const {
      assert(x<size_x);
      assert(y<size_y);
      return array[y*size_x+x];
    }
    
  private:
    T* array;
    size_t size_x, size_y;
  };
}
