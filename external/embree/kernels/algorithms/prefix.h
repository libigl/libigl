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

namespace embree
{
  /*! Implementation of parallel prefix operations (e.g. parallel
   *  prefix sums). The prefix operations start with the identity
   *  (e.g. zero for prefix sums). */
  template<typename SrcArray, typename DstArray, typename Ty, typename Op>
  class ParallelPrefixOp
  {
#if defined(__MIC__)
    static const size_t MAX_THREADS = MAX_MIC_THREADS;
    static const size_t SINGLE_THREAD_THRESHOLD = 50000;
#else
    static const size_t MAX_THREADS = 32;
    static const size_t SINGLE_THREAD_THRESHOLD = 3000000;
#endif

  public:
    ParallelPrefixOp () {}

    class Task
    {
    public:

      Task (ParallelPrefixOp* parent, const SrcArray& src, DstArray& dst, const size_t N, const Op op, const Ty id)
	: parent(parent), src(src), dst(dst), N(N), op(op), id(id)
      {
	/* perform single threaded prefix operation for small N */
	if (N < SINGLE_THREAD_THRESHOLD) {
          size_t sum=0;
	  for (size_t i=0; i<N; sum+=src[i++]) dst[i] = sum;
          parent->value = sum;
	}

	/* perform parallel prefix operation for large N */
	else 
	{
	  LockStepTaskScheduler* scheduler = LockStepTaskScheduler::instance();
	  const size_t numThreads = min(scheduler->getNumThreads(),MAX_THREADS);

	  /* first calculate range for each block */
	  scheduler->dispatchTask(task_sum,this,0,numThreads);
	  
	  /* now calculate prefix_op for each block */
	  scheduler->dispatchTask(task_prefix_op,this,0,numThreads);
	}
      }
      
      /* task that sums up a block of elements */
      void sum(const size_t threadIndex, const size_t threadCount)
      {
	const size_t start = (threadIndex+0)*N/threadCount;
	const size_t end   = (threadIndex+1)*N/threadCount;
	
	Ty sum = id;
	for (size_t i=start; i<end; i++)
	  sum = op(sum,src[i]);
	
	parent->state[threadIndex] = sum;
      }

      /* task that calculates the prefix operation for a block of elements */
      void prefix_op(const size_t threadIndex, const size_t threadCount)
      {
	const size_t start = (threadIndex+0)*N/threadCount;
	const size_t end   = (threadIndex+1)*N/threadCount;
		
	/* calculate start sum for block */
	Ty count = id;
	for (size_t i=0; i<threadIndex; i++)
	  count += parent->state[i];
	
	/* calculate prefix sums for the block */
	for (size_t i=start; i<end; i++) 
	{
	  const Ty v = src[i];
	  dst[i] = count;
	  count = op(count,v);
	}

        if (threadIndex == threadCount-1) 
          parent->value = count;
      }
      
      static void task_sum (void* data, const size_t threadIndex, const size_t threadCount) { 
	((Task*)data)->sum(threadIndex,threadCount);                          
      }
      
      static void task_prefix_op (void* data, const size_t threadIndex, const size_t threadCount) { 
	((Task*)data)->prefix_op(threadIndex,threadCount);                          
      }

    private:
      ParallelPrefixOp* const parent;
      const SrcArray& src;               //!< source array
      DstArray& dst;                     //!< destination array
      const size_t N;                    //!< number of elements in the arrays
      const Op op;                       //!< operation to use for prefix sum
      const Ty id;                       //!< identity of the operation
    };

    void operator() (const SrcArray& src, DstArray& dst, const size_t N, const Op op, const Ty id) {
      Task(this,src,dst,N,op,id);
    }

  public:
    Ty state[MAX_THREADS];
    Ty value;
  };

  template<typename Ty>
  struct my_add { // FIXME: use lambda expressions
    Ty operator()(const Ty& a, const Ty& b) const { return a+b; }
  };

  /*! parallel calculation of prefix sums */
  template<typename SrcArray, typename DstArray>
    __forceinline typename SrcArray::value_type parallel_prefix_sum(const SrcArray& src, DstArray& dst, size_t N) 
  {
    typedef typename SrcArray::value_type Value;
    ParallelPrefixOp<SrcArray,DstArray,Value,my_add<Value> > op; op(src,dst,N,my_add<Value>(),0); return op.value;
  }
}
