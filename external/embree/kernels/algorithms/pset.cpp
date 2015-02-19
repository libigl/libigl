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

#include "pset.h"

namespace embree
{
  struct pset_regression_test : public RegressionTest
  {
    pset_regression_test(const char* name) : name(name) {
      registerRegressionTest(this);
    }
    
    bool operator() ()
    {
      bool passed = true;
      printf("%s::%s ... ",TOSTRING(isa),name);
      fflush(stdout);

      /* create vector with random numbers */
      const size_t N = 10000;
      std::vector<uint32> unsorted(N);
      for (size_t i=0; i<N; i++) unsorted[i] = 2*rand();
      
      /* created set from numbers */
      pset<uint32> sorted;
      sorted.init(unsorted);

      /* check that all elements are in the set */
      for (size_t i=0; i<N; i++) {
	passed &= sorted.lookup(unsorted[i]);
      }

      /* check that these elements are not in the set */
      for (size_t i=0; i<N; i++) {
	passed &= !sorted.lookup(unsorted[i]+1);
      }

      /* output if test passed or not */
      if (passed) printf("[passed]\n");
      else        printf("[failed]\n");
      
      return passed;
    }

    const char* name;
  };

  pset_regression_test pset_regression("pset_regression_test");
}
