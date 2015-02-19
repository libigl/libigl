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
#include "common/buffer.h"
#include "algorithms/sort.h"
#include "algorithms/parallel_for.h"

namespace embree
{
  /* implementation of a set of values with parallel construction */
  template<typename T>
  class pset
  {
  public:

    /*! constructors for the parallel set */
    pset () {}
    pset (const std::vector<T>& in) { init(in); }
    pset (const BufferT<T>    & in) { init(in); }

    /*! initialized the parallel set from a vector */
    void init(const std::vector<T>& in) 
    {
      /* reserve sufficient space for all data */
      vec.resize(in.size());
      temp.resize(in.size());

      /* copy data to temporary array */
      parallel_for( size_t(0), in.size(), size_t(4*4096), [&](const range<size_t>& r) 
      {
	for (size_t i=r.begin(); i<r.end(); i++) 
	  vec[i] = in[i];
      });

      /* sort the data */
      radix_sort<T>(vec.data(),temp.data(),vec.size());
    }


    /*! initialized the parallel set from a user buffer */
    void init(const BufferT<T>& in) 
    {
      /* reserve sufficient space for all data */
      vec.resize(in.size());
      temp.resize(in.size());

      /* copy data to temporary array */
      parallel_for( size_t(0), in.size(), size_t(4*4096), [&](const range<size_t>& r) 
      {
	for (size_t i=r.begin(); i<r.end(); i++) 
	  vec[i] = in[i];
      });

      /* parallel radix sort of the data */
      radix_sort<T>(vec.data(),temp.data(),vec.size());
    }

    /*! tests if some element is in the set */
    __forceinline bool lookup(const T& elt) const {
      return std::binary_search(vec.begin(), vec.end(), elt);
    }

    /*! cleans temporary state required for re-construction */
    void cleanup() {
      temp.clear();
    }

    /*! clears all state */
    void clear() {
      vec.clear();
      temp.clear();
    }

  private:
    std::vector<T> vec;   //!< vector containing sorted elements
    std::vector<T> temp;  //!< temporary vector required during construction only
  };
}
