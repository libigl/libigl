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

#ifndef __EMBREE_RANDOM_H__
#define __EMBREE_RANDOM_H__

#include "math.h"

namespace embree
{
  class Random
  {
  public:

    __forceinline Random(const int seed = 27) {
      setSeed(seed);
    }

    __forceinline void setSeed(const int s)
    {
      const int a = 16807;
      const int m = 2147483647;
      const int q = 127773;
      const int r = 2836;
      int j, k;

      if (s == 0) seed = 1;
      else if (s < 0) seed = -s;
      else seed = s;

      for (j = 32+7; j >= 0; j--) {
        k = seed / q;
        seed = a*(seed - k*q) - r*k;
        if (seed < 0) seed += m;
        if (j < 32) table[j] = seed;
      }
      state = table[0];
    }

    __forceinline int getInt()
    {
      const int a = 16807;
      const int m = 2147483647;
      const int q = 127773;
      const int r = 2836;

      int k = seed / q;
      seed = a*(seed - k*q) - r*k;
      if (seed < 0) seed += m;
      int j = state / (1 + (2147483647-1) / 32);
      state = table[j];
      table[j] = seed;

      return state;
    }

    __forceinline int    getInt   (int limit) { return getInt() % limit; }
    __forceinline float  getFloat (         ) { return min(getInt() / 2147483647.0f, 1.0f - float(ulp)); }
    __forceinline double getDouble(         ) { return min(getInt() / 2147483647.0 , 1.0  - double(ulp)); }

  private:
    int seed;
    int state;
    int table[32];
  };
}

#endif
