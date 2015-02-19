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

#include "prefix.h"

namespace embree
{
  struct prefix_sum_regression_test : public RegressionTest
  {
    prefix_sum_regression_test(const char* name) : name(name) {
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
	/* initialize array with random numbers */
        uint32 sum0 = 0;
	std::vector<uint32> src(N);
	for (size_t i=0; i<N; i++) {
	  sum0 += src[i] = rand();
        }
        
	/* calculate parallel prefix sum */
	std::vector<uint32> dst(N);
	memset(dst.data(),0,N*sizeof(uint32));
	
	double t0 = getSeconds();
	for (size_t i=0; i<M; i++) {
	  uint32 sum1 = parallel_prefix_sum(src,dst,N);
          passed &= (sum0 == sum1);
        }
	double t1 = getSeconds();
	printf("%zu/%3.2fM ",N,1E-6*double(N*M)/(t1-t0));
	
	/* check if prefix sum is correct */
	for (size_t i=0, sum=0; i<N; sum+=src[i++])
	  passed &= ((uint32)sum == dst[i]);
      }
      
      /* output if test passed or not */
      if (passed) printf("[passed]\n");
      else        printf("[failed]\n");

      return passed;
    }

    const char* name;
  };

  prefix_sum_regression_test prefix_sum_regression("prefix_sum_regression");
}
