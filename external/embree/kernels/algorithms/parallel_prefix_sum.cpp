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

#include "parallel_prefix_sum.h"

namespace embree
{
  struct parallel_prefix_sum_regression_test : public RegressionTest
  {
    parallel_prefix_sum_regression_test(const char* name) : name(name) {
      registerRegressionTest(this);
    }
    
    bool operator() ()
    {
      bool passed = true;
      printf("%s::%s ... ",TOSTRING(isa),name);
      fflush(stdout);

      /* create vector with random numbers */
      const size_t N = 100;
      std::vector<size_t> array(N);
      std::vector<atomic_t> prefix_sum(N);
      for (size_t j=0; j<N; j++)
	array[j] = rand() % 10;
  
      /* dry run only counts */
      ParallelPrefixSumState<size_t> state;
      size_t S0 = parallel_prefix_sum( state, size_t(0), size_t(N), size_t(1024), size_t(0), [&](const range<size_t>& r, const size_t sum) 
      {
        size_t s = 0;
	for (size_t i=r.begin(); i<r.end(); i++)
          s += array[i];
	
        return s;
      }, [](size_t v0, size_t v1) { return v0+v1; });
      
      /* final run calculates prefix sum */
      size_t S1 = parallel_prefix_sum( state, size_t(0), size_t(N), size_t(1024), size_t(0), [&](const range<size_t>& r, const size_t sum) 
      {
        size_t s = 0;
	for (size_t i=r.begin(); i<r.end(); i++) {
	  prefix_sum[i] = sum+s;
          s += array[i];
        }
        return s;
      }, [](size_t v0, size_t v1) { return v0+v1; });

      /* check calculated prefix sum */
      size_t sum=0;
      for (size_t i=0; i<N; sum+=array[i++]) {
        passed &= (prefix_sum[i] == sum);
      }

      passed &= (S0 == sum);
      passed &= (S1 == sum);

      /* output if test passed or not */
      if (passed) printf("[passed]\n");
      else        printf("[failed]\n");
      
      return passed;
    }

    const char* name;
  };

  parallel_prefix_sum_regression_test parallel_prefix_sum_regression("parallel_prefix_sum_regression_test");
}
