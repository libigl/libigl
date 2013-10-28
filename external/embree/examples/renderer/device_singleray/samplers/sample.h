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

#ifndef __EMBREE_SAMPLE_H__
#define __EMBREE_SAMPLE_H__

#include "../default.h"

namespace embree
{
  /*! Sample representation. A sample has a location and a probability
   *  density it got sampled with. */
  template<typename T> struct Sample
  {
    __forceinline Sample           () {}
    __forceinline Sample           ( const Sample& other ) { value = other.value; pdf = other.pdf; }
    __forceinline Sample& operator=( const Sample& other ) { value = other.value; pdf = other.pdf; return *this; }

    __forceinline Sample (const T& value, float pdf = one) : value(value), pdf(pdf) {}

    __forceinline operator const T&( ) const { return value; }
    __forceinline operator       T&( )       { return value; }

  public:
    T value;    //!< location of the sample
    float pdf;  //!< probability density of the sample
  };

  /*! output operator */
  template<typename T> __forceinline std::ostream& operator<<(std::ostream& cout, const Sample<T>& s) {
    return cout << "{ value = " << s.value << ", pdf = " << s.pdf << "}";
  }

  /*! shortcuts for common sample types */
  typedef Sample<int>   Sample1i;
  typedef Sample<float> Sample1f;
  typedef Sample<Vec2f> Sample2f;
  typedef Sample<Vector3f> Sample3f;
}

#endif
