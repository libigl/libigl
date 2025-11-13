// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PARALLEL_FOR_H
#define IGL_PARALLEL_FOR_H
#include "igl_inline.h"
#include <functional>

//#warning "Defining IGL_PARALLEL_FOR_FORCE_SERIAL"
//#define IGL_PARALLEL_FOR_FORCE_SERIAL

namespace igl
{
  /// Functional implementation of a basic, open-mp style, parallel
  /// for loop. If the inner block of a for-loop can be rewritten/encapsulated in
  /// a single (anonymous/lambda) function call `func` so that the serial code
  /// looks like:
  ///
  /// \code{cpp}
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       func(i);
  ///     }
  /// \endcode
  ///
  /// then `parallel_for(loop_size,func,min_parallel)` will use as many threads as
  /// available on the current hardware to parallelize this for loop so long as
  /// loop_size<min_parallel, otherwise it will just use a serial for loop.
  ///
  /// Often if your code looks like:
  ///
  /// \code{cpp}
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       …
  ///     }
  /// \endcode
  ///
  /// Then you can make a minimal two-line change to parallelize it:
  ///
  /// \code{cpp}
  ///     //for(int i = 0;i<loop_size;i++)
  ///     parallel_for(loop_size,[&](int i)
  ///     {
  ///       …
  ///     }
  ///     ,1000);
  /// \endcode
  ///
  /// @param[in] loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
  /// @param[in] func  function handle taking iteration index as only argument to compute
  ///     inner block of for loop I.e. for(int i ...){ func(i); }
  /// @param[in] min_parallel  min size of loop_size such that parallel (non-serial)
  ///     thread pooling should be attempted {0}
  /// @return true iff thread pool was invoked
  template<typename Index, typename FunctionType >
  inline bool parallel_for(
    const Index loop_size,
    const FunctionType & func,
    const size_t min_parallel=0);

  /// Functional implementation of an open-mp style, parallel for loop with
  /// accumulation. For example, serial code separated into n chunks (each to be
  /// parallelized with a thread) might look like:
  ///
  /// \code{cpp}
  ///     Eigen::VectorXd S;
  ///     const auto & prep_func = [&S](int n){ S = Eigen:VectorXd::Zero(n); };
  ///     const auto & func = [&X,&S](int i, int t){ S(t) += X(i); };
  ///     const auto & accum_func = [&S,&sum](int t){ sum += S(t); };
  ///     prep_func(n);
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       func(i,i%n);
  ///     }
  ///     double sum = 0;
  ///     for(int t = 0;t<n;t++)
  ///     {
  ///       accum_func(t);
  ///     }
  /// \endcode
  ///
  /// @param[in] loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
  /// @param[in] prep_func function handle taking n >= number of threads as only
  ///     argument
  /// @param[in] func  function handle taking iteration index i and thread id t as only
  ///     arguments to compute inner block of for loop I.e.
  ///     for(int i ...){ func(i,t); }
  /// @param[in] accum_func  function handle taking thread index as only argument, to be
  ///     called after all calls of func, e.g., for serial accumulation across
  ///     all n (potential) threads, see n in description of prep_func.
  /// @param[in] min_parallel  min size of loop_size such that parallel (non-serial)
  ///     thread pooling should be attempted {0}
  /// @return true iff thread pool was invoked
  template<
    typename Index,
    typename PrepFunctionType,
    typename FunctionType,
    typename AccumFunctionType
    >
  inline bool parallel_for(
    const Index loop_size,
    const PrepFunctionType & prep_func,
    const FunctionType & func,
    const AccumFunctionType & accum_func,
    const size_t min_parallel=0);
}

// Implementation

#include "default_num_threads.h"

#include <cmath>
#include <cassert>
#include <thread>
#include <vector>
#include <algorithm>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace igl {
namespace internal
{

inline bool & worker_flag()
{
  static thread_local bool flag = false;
  return flag;
}

inline bool is_worker_thread()
{
  return worker_flag();
}

inline void set_worker_thread(bool v)
{
  worker_flag() = v;
}

} // namespace internal
} // namespace igl

namespace igl {
namespace internal
{

// Simple shared thread pool using only std::
class ThreadPool
{
public:
  using Task = std::function<void()>;

  // First call fixes the pool size; later calls ignore nthreads_hint.
  static ThreadPool & instance(size_t nthreads_hint)
  {
    static ThreadPool pool(nthreads_hint);
    return pool;
  }

  size_t size() const
  {
    return workers.size();
  }

  void enqueue(Task task)
  {
    {
      std::unique_lock<std::mutex> lock(mutex);
      tasks.emplace(std::move(task));
    }
    cv.notify_one();
  }

private:
  ThreadPool(size_t nthreads_hint)
    : stop(false)
  {
    size_t nthreads = nthreads_hint == 0 ? 1 : nthreads_hint;
    workers.reserve(nthreads);
    for(size_t i = 0; i < nthreads; ++i)
    {
      workers.emplace_back([this]()
          {
          igl::internal::set_worker_thread(true); // <- mark this thread as a pool worker

          for(;;)
          {
          Task task;
          {
          std::unique_lock<std::mutex> lock(mutex);
          cv.wait(lock, [this]()
              {
              return stop || !tasks.empty();
              });

          if(stop && tasks.empty())
          {
          return;
          }

          task = std::move(tasks.front());
          tasks.pop();
          }
          task();
          }
          });
    }
  }

  ~ThreadPool()
  {
    {
      std::unique_lock<std::mutex> lock(mutex);
      stop = true;
    }
    cv.notify_all();
    for(std::thread & t : workers)
    {
      if(t.joinable()) t.join();
    }
  }

  std::vector<std::thread> workers;
  std::queue<Task> tasks;
  mutable std::mutex mutex;
  std::condition_variable cv;
  bool stop;
};

} // namespace internal
} // namespace igl

template<typename Index, typename FunctionType >
inline bool igl::parallel_for(
  const Index loop_size,
  const FunctionType & func,
  const size_t min_parallel)
{
  // no-op preparation/accumulation
  const auto & no_op = [](const size_t /*n_or_t*/){};
  // two-parameter wrapper ignoring thread id
  const auto & wrapper = [&func](Index i, size_t /*t*/){ func(i); };
  return parallel_for(loop_size, no_op, wrapper, no_op, min_parallel);
}


template<
  typename Index,
  typename PreFunctionType,
  typename FunctionType,
  typename AccumFunctionType>
inline bool igl::parallel_for(
  const Index loop_size,
  const PreFunctionType & prep_func,
  const FunctionType & func,
  const AccumFunctionType & accum_func,
  const size_t min_parallel)
{
  assert(loop_size >= 0);
  if (loop_size == 0) return false;

  // If we're already inside a ThreadPool worker, run serial to avoid nested
  // deadlock with the global pool.
  if (igl::internal::is_worker_thread())
  {
    prep_func(1);
    for (Index i = 0; i < loop_size; ++i)
    {
      func(i, 0);
    }
    accum_func(0);
    return false;
  }

#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
  const size_t configured_threads = 1;
#else
  const size_t configured_threads = igl::default_num_threads();
#endif

  if (loop_size < static_cast<Index>(min_parallel) || configured_threads <= 1)
  {
    // Serial fallback
    prep_func(1);
    for (Index i = 0; i < loop_size; ++i)
    {
      func(i, 0);
    }
    accum_func(0);
    return false;
  }

  // --- Parallel branch using shared thread pool ---

  auto & pool = igl::internal::ThreadPool::instance(configured_threads);
  const size_t pool_threads = std::max<size_t>(1, pool.size());

  // Match old semantics: prep called with number of *potential* threads.
  prep_func(pool_threads);

  // Number of "logical jobs" (chunks of the index range).
  const size_t jobs = static_cast<size_t>(
    std::min<Index>(loop_size, static_cast<Index>(pool_threads)));

  struct Group
  {
    std::mutex mutex;
    std::condition_variable cv;
    std::atomic<size_t> remaining;
  };

  auto group = std::make_shared<Group>();
  group->remaining.store(jobs, std::memory_order_relaxed);

  const Index total = loop_size;
  const Index base  = total / static_cast<Index>(jobs);
  const Index rem   = total % static_cast<Index>(jobs);

  for (size_t t = 0; t < jobs; ++t)
  {
    const Index start =
      static_cast<Index>(t) * base
      + std::min<Index>(static_cast<Index>(t), rem);

    const Index end = start + base + (t < static_cast<size_t>(rem) ? 1 : 0);

    pool.enqueue([group, &func, start, end, t]()
    {
      // Each job processes its contiguous slice [start, end)
      for (Index k = start; k < end; ++k)
      {
        func(k, t);
      }

      // Signal completion of this job.
      if (group->remaining.fetch_sub(1, std::memory_order_acq_rel) == 1)
      {
        std::unique_lock<std::mutex> lock(group->mutex);
        group->cv.notify_one();
      }
    });
  }

  // Wait for all jobs for this parallel_for call to finish.
  {
    std::unique_lock<std::mutex> lock(group->mutex);
    group->cv.wait(lock, [&group]()
    {
      return group->remaining.load(std::memory_order_acquire) == 0;
    });
  }

  // Accumulate across all potential threads (same as original implementation).
  for (size_t t = 0; t < pool_threads; ++t)
  {
    accum_func(t);
  }

  return true;
}

#endif

