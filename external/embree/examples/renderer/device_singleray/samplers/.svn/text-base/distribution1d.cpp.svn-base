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

#include "distribution1d.h"
#include <algorithm>

namespace embree
{
  Distribution1D::Distribution1D()
    : size(0), PDF(NULL), CDF(NULL) {}

  Distribution1D::Distribution1D(const float* f, const size_t size_in) {
    init(f, size_in);
  }

  Distribution1D::~Distribution1D() {
    if (PDF) delete[] PDF; PDF = NULL;
    if (CDF) delete[] CDF; CDF = NULL;
  }

  void Distribution1D::readback(float *cdf, float *pdf)
  {
    assert(CDF);
    assert(PDF);
    std::copy(CDF,CDF+size+1,cdf);
    std::copy(PDF,PDF+size,pdf);
  }

  void Distribution1D::init(const float* f, const size_t size_in)
  {
    /*! create arrays */
    size = size_in;
    PDF = new float[size];
    CDF = new float[size+1];

    /*! accumulate the function f */
    CDF[0] = 0.0f;
    for (size_t i=1; i<size+1; i++)
      CDF[i] = CDF[i-1] + f[i-1];

    /*! compute reciprocal sum */
    float rcpSum = CDF[size] == 0.0f ? 0.0f : rcp(CDF[size]);

    /*! normalize the probability distribution and cumulative distribution */
    for (size_t i = 1; i<size+1; i++) {
      PDF[i-1] = f[i-1] * rcpSum * size;
      CDF[i] *= rcpSum;
     }
    CDF[size] = 1.0f;
  }

  Sample1f Distribution1D::sample(const float u) const
  {
    /*! coarse sampling of the distribution */
    float* pointer = std::upper_bound(CDF, CDF+size, u);
    int index = clamp(int(pointer-CDF-1),0,int(size)-1);
    
    /*! refine sampling linearly by assuming the distribution being a step function */
    float fraction = (u - CDF[index]) * rcp(CDF[index+1] - CDF[index]);
    return Sample1f(float(index)+fraction,PDF[index]);
  }

  float Distribution1D::pdf(const float p) const {
    return PDF[clamp(int(p*size),0,int(size)-1)];
  }
}

