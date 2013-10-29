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

#ifndef __EMBREE_PERMUTATION_H__
#define __EMBREE_PERMUTATION_H__

#include "sys/platform.h"
#include "random.h"

#include <algorithm>

namespace embree
{
  /*! Random permutation. */
  class Permutation
  {
  public:

    /*! Creates a random permutation. Uses system random number generator. */
    Permutation(int size)
    {
      elts = size;
      perm = new int[elts];
      for (int i=0; i<elts; i++) perm[i] = i;
      for (int i=0; i<elts; i++) std::swap(perm[i],perm[random<uint32>()%elts]);
    }

    /*! Creates a random permutation. Random number generator is provided as argument. */
    Permutation(int size, Random& rng)
    {
      elts = size;
      perm = new int[elts];
      for (int i=0; i<elts; i++) perm[i] = i;
      for (int i=0; i<elts; i++) std::swap(perm[i],perm[rng.getInt(elts)]);
    }

    /*! Destroys the permutation. */
    ~Permutation() {
      if (perm) delete[] perm; perm = NULL;
    }

    /*! Returns the size of the permutation. */
    __forceinline int size() const {
      return elts;
    }

    /*! Returns the i'th element of the permutation. */
    __forceinline int operator[](int i) const {
      assert(i >= 0 && i < elts);
      return perm[i];
    }

  private:
    int elts;    //!< Size of the permutation.
    int* perm;   //!< Array storing the permutation.
  };
}

#endif
