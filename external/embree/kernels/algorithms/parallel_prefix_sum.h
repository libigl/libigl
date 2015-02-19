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

#include "parallel_for.h"

namespace embree
{
  template<typename Value>
    struct ParallelPrefixSumState 
  {
    enum { MAX_TASKS = 32 };
    Value counts[MAX_TASKS];
    Value sums  [MAX_TASKS];
  };

  template<typename Index, typename Value, typename Func, typename Reduction>
    __forceinline Value parallel_prefix_sum( ParallelPrefixSumState<Value>& state, Index first, Index last, Index minStepSize, const Value& identity, const Func& func, const Reduction& reduction)
  {
    /* calculate number of tasks to use */
    LockStepTaskScheduler* scheduler = LockStepTaskScheduler::instance();
    const size_t numThreads = scheduler->getNumThreads();
    const size_t numBlocks  = (last-first+minStepSize-1)/minStepSize;
    const size_t taskCount  = min(numThreads,numBlocks,size_t(ParallelPrefixSumState<Value>::MAX_TASKS));

    /* perform parallel prefix sum */
    parallel_for(taskCount, [&](const size_t taskIndex)
    {
      const size_t i0 = first+(taskIndex+0)*(last-first)/taskCount;
      const size_t i1 = first+(taskIndex+1)*(last-first)/taskCount;
      state.counts[taskIndex] = func(range<size_t>(i0,i1),state.sums[taskIndex]);
    });

    /* calculate prefix sum */
    Value sum=0;
    for (size_t i=0; i<taskCount; i++)
    {
      const Value c = state.counts[i];
      state.sums[i] = sum;
      sum=reduction(sum,c);
    }

    return sum;
  }
}
