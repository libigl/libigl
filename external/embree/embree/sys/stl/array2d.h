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

#ifndef __EMBREE_ARRAY2D_H__
#define __EMBREE_ARRAY2D_H__

#include "../platform.h"

namespace embree
{
  template<typename T>
  class Array2D
  {
  public:
    Array2D() : sizeX(0), sizeY(0), data(NULL) {}

    Array2D(size_t sizeX, size_t sizeY) : sizeX(sizeX), sizeY(sizeY) {
      data = new T*[sizeX];
      for (size_t x = 0; x < sizeX; x++)
        data[x] = new T[sizeY];
    }

    ~Array2D() {
      for (size_t x = 0; x < sizeX; x++)
        delete[] data[x];
      delete[] data;
    }

    operator const T**() const { return const_cast<const T**>(data); }
    operator T**() { return data; }
    const T& get(const size_t x, const size_t y) const { return data[x][y]; }
    T& get(const size_t x, const size_t y) { return data[x][y]; }
    void set(const size_t x, const size_t y, const T& value) { data[x][y] = value; }

  private:
    size_t sizeX;
    size_t sizeY;
    T** data;
  };
}

#endif
