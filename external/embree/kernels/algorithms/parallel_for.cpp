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

#include "parallel_for.h"

namespace embree
{
  struct parallel_for_regression_test : public RegressionTest
  {
    parallel_for_regression_test(const char* name) : name(name) {
      registerRegressionTest(this);
    }
    
    bool operator() ()
    {
      bool passed = true;
      printf("%s::%s ... ",TOSTRING(isa),name);
      fflush(stdout);

      const size_t M = 10;
      for (size_t N=10; N<10000000; N*=2.1f)
      {
        /* sequentially calculate sum of squares */
        size_t sum0 = 0;
        for (size_t i=0; i<N; i++) {
          sum0 += i*i;
        }

        /* parallel calculation of sum of squares */
	double t0 = getSeconds();
        for (size_t m=0; m<M; m++)
        {
          AtomicCounter sum1 = 0;
          parallel_for( size_t(0), size_t(N), size_t(1024), [&](const range<size_t>& r) 
          {
            size_t s = 0;
            for (size_t i=r.begin(); i<r.end(); i++) 
              s += i*i;
            sum1 += s;
          });
          passed = sum0 == sum1;
        }
	double t1 = getSeconds();
	printf("%zu/%3.2fM ",N,1E-6*double(N*M)/(t1-t0));
      }
      
      /* output if test passed or not */
      if (passed) printf("[passed]\n");
      else        printf("[failed]\n");
      
      return passed;
    }

    const char* name;
  };

  parallel_for_regression_test parallel_for_regression("parallel_for_regression_test");
}
