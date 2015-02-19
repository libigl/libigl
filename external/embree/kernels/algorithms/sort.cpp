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

#include "sort.h"

namespace embree
{
  //ParallelRadixSort shared_radix_sort_state;
  //ParallelRadixSortT<uint32> radix_sort_u32(shared_radix_sort_state);
  //ParallelRadixSortT<uint64> radix_sort_u64(shared_radix_sort_state);
  
  template<typename Key>
  struct RadixSortRegressionTest : public RegressionTest
  {
    RadixSortRegressionTest(const char* name) : name(name) {
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
	std::vector<Key> src(N); memset(src.data(),0,N*sizeof(Key));
	std::vector<Key> tmp(N); memset(tmp.data(),0,N*sizeof(Key));
	for (size_t i=0; i<N; i++) src[i] = uint64(rand())*uint64(rand());
	
	/* calculate checksum */
	Key sum0 = 0; for (size_t i=0; i<N; i++) sum0 += src[i];
        
	/* sort numbers */
	double t0 = getSeconds();
	for (size_t i=0; i<M; i++) {
          radix_sort<Key>(src.data(),tmp.data(),N);
        }
	double t1 = getSeconds();
	printf("%zu/%3.2fM ",N,1E-6*double(N*M)/(t1-t0));
	
	/* calculate checksum */
	Key sum1 = 0; for (size_t i=0; i<N; i++) sum1 += src[i];
	if (sum0 != sum1) passed = false;
        
	/* check if numbers are sorted */
	for (size_t i=1; i<N; i++)
	  passed &= src[i-1] <= src[i];
      }
      
      /* output if test passed or not */
      if (passed) printf("[passed]\n");
      else        printf("[failed]\n");

      return passed;
    }

    const char* name;
  };

  RadixSortRegressionTest<uint32> test_u32("RadixSortRegressionTestU32");
  RadixSortRegressionTest<uint64> test_u64("RadixSortRegressionTestU64");
}
