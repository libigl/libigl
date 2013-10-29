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

#ifndef __EMBREE_PATTERNS_H__
#define __EMBREE_PATTERNS_H__

#include "default.h"
#include "sys/stl/array2d.h"
#include "math/permutation.h"

namespace embree
{
  /*! Create a set of n jittered 1D samples, using the provided 
   *  random number generator. */
  __forceinline void jittered(float* samples, const uint32 n, Random& rng)
  {
    float scale = 1.0f / n;
    Permutation permutation(n, rng);
    for (uint32 i = 0; i < n; i++) {
      samples[permutation[i]] = (float(i) + rng.getFloat()) * scale;
    }
  }

  /*! Create a set of n multi-jittered 2D samples, using the provided 
   *  random number generator. */
  void multiJittered(Vec2f* samples, const uint32 N, Random& rng)
  {
    uint32 b = (uint32)sqrtf(float(N));
    if (b*b < N) b++;
    Array2D<Vec2f> grid(b,b);

    vector_t<uint32> numbers(b);
    for (uint32 i = 0; i < b; i++)
      numbers[i] = i;

    for (uint32 i = 0; i < b; i++) {
      numbers.shuffle(rng);
      for (uint32 j = 0; j < b; j++) {
        ((Vec2f**)grid)[i][j][0] = float(i)/float(b) + (numbers[j]+rng.getFloat())/float(b*b);
      }
    }

    for (uint32 i = 0; i < b; i++) {
      numbers.shuffle(rng);
      for (uint32 j = 0; j < b; j++) {
        ((Vec2f**)grid)[j][i][1] = float(i)/float(b) + (numbers[j]+rng.getFloat())/float(b*b);
      }
    }

    Permutation permutation(N, rng);
    for (uint32 n=0; n<N; n++) {
      uint32 np = permutation[n];
      samples[n] = grid.get(np/b,np%b);
    }
  }
}

#endif
