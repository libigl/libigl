#include <test_common.h>
#include <igl/parallel_for.h>
#include <atomic>
#include <vector>
#include <numeric>
#include <thread>
#include <set>

TEST_CASE("parallel_for: serial_fallback", "[igl][parallel_for]")
{
  // loop_size < min_parallel ⇒ must run serial
  std::vector<int> vals(10, 0);

  bool used_parallel = igl::parallel_for(
    (int)vals.size(),
    [&](int i){ vals[i] = 1; },
    /*min_parallel=*/1000
  );

  REQUIRE(used_parallel == false);
  for (int v : vals)
    REQUIRE(v == 1);
}

TEST_CASE("parallel_for: basic_parallelism", "[igl][parallel_for]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const int N = 20000;
  std::vector<int> hit(N, 0);
  std::atomic<int> counter(0);

  bool used_parallel = igl::parallel_for(
    N,
    [&](int i)
    {
      hit[i] = 1;
      counter.fetch_add(1, std::memory_order_relaxed);
    },
    /*min_parallel=*/1
  );

  REQUIRE(used_parallel == true);
  REQUIRE(counter.load() == N);
  for (int v : hit)
    REQUIRE(v == 1);
}

TEST_CASE("parallel_for: accumulation", "[igl][parallel_for]")
{
  const int N = 10000;

  // Per-thread buckets
  std::vector<double> buckets;

  const auto prep = [&](size_t nt)
  {
    buckets.assign(nt, 0.0);
  };

  const auto func = [&](int /*i*/, size_t t)
  {
    buckets[t] += 1.0; // increment per-thread
  };

  double total = 0.0;
  const auto accum = [&](size_t t)
  {
    total += buckets[t];
  };

  bool used_parallel = igl::parallel_for(N, prep, func, accum, 1);

  if(igl::default_num_threads() > 1) 
  { 
    REQUIRE(used_parallel == true);
  }
  REQUIRE(total == Approx((double)N));
}

TEST_CASE("parallel_for: equivalence_to_serial", "[igl][parallel_for]")
{
  const int N = 15000;

  // serial result
  std::vector<int> S(N);
  for (int i = 0; i < N; ++i) S[i] = i*i;

  // parallel result
  std::vector<int> P(N);

  igl::parallel_for(
    N,
    [&](int i){ P[i] = i*i; },
    /*min_parallel=*/1
  );

  for (int i = 0; i < N; ++i)
    REQUIRE(P[i] == S[i]);
}

TEST_CASE("parallel_for: min_parallel_threshold", "[igl][parallel_for]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const int N = 500;
  std::vector<int> A(N,0), B(N,0);

  bool p1 = igl::parallel_for(
    N, [&](int i){ A[i] = i; },
    /*min_parallel=*/1000 // too large → serial
  );
  bool p2 = igl::parallel_for(
    N, [&](int i){ B[i] = i; },
    /*min_parallel=*/1 // small → parallel
  );

  REQUIRE(p1 == false);
  REQUIRE(p2 == true);
  REQUIRE(A == B);
}

TEST_CASE("parallel_for: nested_calls", "[igl][parallel_for]")
{
  const int N = 2000;
  std::vector<int> out(N, 0);

  bool used_parallel = igl::parallel_for(
    N,
    [&](int i)
    {
      // a tiny nested parallel_for
      igl::parallel_for(
        10,
        [&](int j){ (void)j; }, 1
      );
      out[i] = 1;
    },
    /*min_parallel=*/1
  );

  if(igl::default_num_threads() > 1) 
  {
    REQUIRE(used_parallel == true);
  }
  for (int v : out)
    REQUIRE(v == 1);
}

// -----------------------------------------------------------------------------
// Additional tests
// -----------------------------------------------------------------------------

TEST_CASE("parallel_for: zero_iterations_does_nothing", "[igl][parallel_for]")
{
  std::atomic<int> prep_calls(0);
  std::atomic<int> func_calls(0);
  std::atomic<int> accum_calls(0);

  const auto prep = [&](size_t /*nt*/)
  {
    prep_calls.fetch_add(1, std::memory_order_relaxed);
  };

  const auto func = [&](int /*i*/, size_t /*t*/)
  {
    func_calls.fetch_add(1, std::memory_order_relaxed);
  };

  const auto accum = [&](size_t /*t*/)
  {
    accum_calls.fetch_add(1, std::memory_order_relaxed);
  };

  bool used_parallel = igl::parallel_for(
    0,
    prep,
    func,
    accum,
    /*min_parallel=*/1
  );

  REQUIRE(used_parallel == false);
  REQUIRE(prep_calls.load() == 0);
  REQUIRE(func_calls.load() == 0);
  REQUIRE(accum_calls.load() == 0);
}

TEST_CASE("parallel_for: min_parallel_equal_threshold", "[igl][parallel_for]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const int N = 1024;
  std::vector<int> A(N,0), B(N,0);

  // spec: serial if loop_size < min_parallel, so == should be parallel
  bool serial_used = igl::parallel_for(
    N, [&](int i){ A[i] = i; },
    /*min_parallel=*/N+1
  );
  bool parallel_used = igl::parallel_for(
    N, [&](int i){ B[i] = i; },
    /*min_parallel=*/N
  );

  REQUIRE(serial_used == false);
  REQUIRE(parallel_used == true);
  REQUIRE(A == B);
}

TEST_CASE("parallel_for: thread_id_range_and_accum_calls", "[igl][parallel_for]")
{
  const int N = 10000;

  std::vector<long long> bucket_counts;
  std::atomic<size_t> prep_nt(0);
  std::atomic<size_t> max_t_seen(0);
  std::atomic<size_t> accum_calls(0);

  const auto prep = [&](size_t nt)
  {
    prep_nt.store(nt, std::memory_order_relaxed);
    bucket_counts.assign(nt, 0); // plain ints, no atomics needed
  };

  const auto func = [&](int /*i*/, size_t t)
  {
    // track max t across threads
    size_t cur = max_t_seen.load(std::memory_order_relaxed);
    while (t > cur && !max_t_seen.compare_exchange_weak(
                        cur, t,
                        std::memory_order_relaxed,
                        std::memory_order_relaxed))
    {
      // spin until we either win or see a newer/bigger cur
    }

    // Each logical t corresponds to a single job in this implementation,
    // so all calls with the same t are executed on the same worker thread,
    // sequentially. No race here; we can use a plain increment.
    bucket_counts[t] += 1;
  };

  const auto accum = [&](size_t /*t*/)
  {
    accum_calls.fetch_add(1, std::memory_order_relaxed);
  };

  bool used_parallel = igl::parallel_for(
    N, prep, func, accum,
    /*min_parallel=*/1
  );


  const size_t nt = prep_nt.load();
  REQUIRE(nt >= 1);

  // t must always be < nt
  REQUIRE(max_t_seen.load() < nt);

  // accum must be called once per potential thread
  REQUIRE(accum_calls.load() == nt);

  // Sanity: total counted iterations == N
  long long total = 0;
  for (size_t t = 0; t < nt; ++t)
    total += bucket_counts[t];
  REQUIRE(total == N);
}



TEST_CASE("parallel_for: nested_inner_serial_fallback", "[igl][parallel_for]")
{
  const int N = 1000;
  std::vector<int> outer_hits(N, 0);

  std::atomic<bool> inner_parallel_seen(false);

  bool outer_parallel = igl::parallel_for(
    N,
    [&](int i)
    {
      bool inner_used_parallel = igl::parallel_for(
        10,
        [&](int j){ (void)j; },
        /*min_parallel=*/1
      );

      if (inner_used_parallel)
      {
        inner_parallel_seen.store(true, std::memory_order_relaxed);
      }

      outer_hits[i] = 1;
    },
    /*min_parallel=*/1
  );

  if(igl::default_num_threads() > 1) 
  { 
    REQUIRE(outer_parallel == true);
  }
  for (int v : outer_hits)
    REQUIRE(v == 1);

  // With the is_worker_thread() guard in the implementation,
  // inner calls from pool workers should always be serial.
  REQUIRE(inner_parallel_seen.load() == false);
}



TEST_CASE("parallel_for: deep_nested_calls", "[igl][parallel_for]")
{
  const int N = 256;
  std::vector<int> hits(N, 0);

  bool outer_parallel = igl::parallel_for(
    N,
    [&](int i)
    {
      igl::parallel_for(
        8,
        [&](int j)
        {
          // third level
          igl::parallel_for(
            4,
            [&](int k){ (void)k; },
            1
          );
          (void)j;
        },
        1
      );
      hits[i] = 1;
    },
    /*min_parallel=*/1
  );

  if(igl::default_num_threads() > 1) 
  { 
    REQUIRE(outer_parallel == true);
  }
  for (int v : hits)
    REQUIRE(v == 1);
}

TEST_CASE("parallel_for: many_small_jobs_reuse_pool", "[igl][parallel_for]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const int iterations = 200;
  const int N = 64;

  std::vector<int> buf(N);

  for (int it = 0; it < iterations; ++it)
  {
    std::fill(buf.begin(), buf.end(), 0);

    bool used_parallel = igl::parallel_for(
      N,
      [&](int i){ buf[i] = it; },
      /*min_parallel=*/1
    );
    if(igl::default_num_threads() > 1) 
    { 
      REQUIRE(used_parallel == true);
    }

    for (int i = 0; i < N; ++i)
      REQUIRE(buf[i] == it);
  }
}

TEST_CASE("parallel_for: different_index_types", "[igl][parallel_for]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const long long N = 12345;

  std::vector<int> buf((size_t)N, 0);

  bool used_parallel = igl::parallel_for(
    N,
    [&](long long i)
    {
      buf[(size_t)i] = 1;
    },
    /*min_parallel=*/1
  );

  REQUIRE(used_parallel == true);
  for (int v : buf)
    REQUIRE(v == 1);
}

TEST_CASE("parallel_for: accumulation_equivalence_to_serial_sum", "[igl][parallel_for]")
{
  const int N = 10000;

  // serial sum
  long long serial_sum = 0;
  for (int i = 0; i < N; ++i)
  {
    serial_sum += i;
  }

  // parallel sum: S[t] accumulates partial sums, then accum collects.
  std::vector<long long> S;
  long long parallel_sum = 0;

  const auto prep = [&](size_t nt)
  {
    S.assign(nt, 0);
  };

  const auto func = [&](int i, size_t t)
  {
    S[t] += i;
  };

  const auto accum = [&](size_t t)
  {
    parallel_sum += S[t];
  };

  bool used_parallel = igl::parallel_for(
    N, prep, func, accum,
    /*min_parallel=*/1
  );

  if(igl::default_num_threads() > 1)
  {
    REQUIRE(used_parallel == true);
  }
  REQUIRE(parallel_sum == serial_sum);
}

#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
TEST_CASE("parallel_for: force_serial_macro", "[igl][parallel_for]")
{
  // If compiled with IGL_PARALLEL_FOR_FORCE_SERIAL, we must never see parallel.
  const int N = 1000;
  std::vector<int> buf(N, 0);

  bool used_parallel = igl::parallel_for(
    N,
    [&](int i){ buf[i] = i; },
    /*min_parallel=*/1
  );

  REQUIRE(used_parallel == false);
  for (int i = 0; i < N; ++i)
    REQUIRE(buf[i] == i);
}
#endif

//#define IGL_PARALLEL_FOR_TIMING_TESTS
#ifdef IGL_PARALLEL_FOR_TIMING_TESTS

#include <chrono>
#include <cmath>

// Helper alias
using igl_pf_clock = std::chrono::steady_clock;

TEST_CASE("parallel_for: timing_large_loop", "[igl][parallel_for][timing]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  const int N = 5'000'000;

  std::vector<double> a(N), b(N);
  for (int i = 0; i < N; ++i)
  {
    a[i] = 0.5 * i;
  }

  // --- Serial baseline ---
  auto serial_start = igl_pf_clock::now();
  for (int i = 0; i < N; ++i)
  {
    // mildly non-trivial work to avoid being optimized away
    b[i] = std::sqrt(a[i] * a[i] + 1.0);
  }
  auto serial_end = igl_pf_clock::now();

  auto serial_ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(
      serial_end - serial_start).count();

  // --- Parallel version ---
  std::fill(b.begin(), b.end(), 0.0);

  auto parallel_start = igl_pf_clock::now();
  bool used_parallel = igl::parallel_for(
    N,
    [&](int i)
    {
      b[i] = std::sqrt(a[i] * a[i] + 1.0);
    },
    /*min_parallel=*/1
  );
  auto parallel_end = igl_pf_clock::now();

  auto parallel_ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(
      parallel_end - parallel_start).count();

  INFO("timing_large_loop: serial_ms   = " << serial_ms);
  INFO("timing_large_loop: parallel_ms = " << parallel_ms);
  INFO("timing_large_loop: used_parallel = " << used_parallel);

  // Sanity: results should match a re-run of serial
  std::vector<double> c(N);
  for (int i = 0; i < N; ++i)
  {
    c[i] = std::sqrt(a[i] * a[i] + 1.0);
  }
  for (int i = 0; i < N; ++i)
  {
    REQUIRE(b[i] == Approx(c[i]));
  }

  // Very soft performance assertion:
  // If we actually ran in parallel and the serial baseline took > 0 ms,
  // then parallel should not be crazy slower (e.g., 10x).
  if (used_parallel && serial_ms > 0)
  {
    double ratio = (parallel_ms > 0)
      ? double(parallel_ms) / double(serial_ms)
      : 0.0;

    INFO("timing_large_loop: parallel / serial ratio = " << ratio);

    // Soft bound: allow parallel to be up to 10x slower in worst case.
    CHECK(ratio < 10.0);
  }
}

TEST_CASE("parallel_for: timing_many_small_jobs", "[igl][parallel_for][timing]")
{
  if(igl::default_num_threads() <= 1) { SUCCEED("Only one hardware thread; nested parallel test skipped."); return; }
  // This is meant to stress the thread pool reuse behavior: many small jobs.
  const int iterations = 500;
  const int N = 1024;

  std::vector<double> data(N, 1.0);

  // --- Serial: do all work in a single loop ---
  auto serial_start = igl_pf_clock::now();
  double serial_sum = 0.0;
  for (int it = 0; it < iterations; ++it)
  {
    for (int i = 0; i < N; ++i)
    {
      serial_sum += data[i] * 0.5;
    }
  }
  auto serial_end = igl_pf_clock::now();

  auto serial_ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(
      serial_end - serial_start).count();

  // --- Parallel: same total work, but split into many parallel_for calls ---
  auto parallel_start = igl_pf_clock::now();
  double parallel_sum = 0.0;

  for (int it = 0; it < iterations; ++it)
  {
    double local_sum = 0.0;

    // Here we use the accum-variant to test that path too.
    std::vector<double> buckets;
    const auto prep = [&](size_t nt)
    {
      buckets.assign(nt, 0.0);
    };
    const auto func = [&](int i, size_t t)
    {
      buckets[t] += data[i] * 0.5;
    };
    const auto accum = [&](size_t t)
    {
      local_sum += buckets[t];
    };

    (void)igl::parallel_for(
      N, prep, func, accum,
      /*min_parallel=*/1
    );

    parallel_sum += local_sum;
  }

  auto parallel_end = igl_pf_clock::now();

  auto parallel_ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(
      parallel_end - parallel_start).count();

  INFO("timing_many_small_jobs: serial_ms   = " << serial_ms);
  INFO("timing_many_small_jobs: parallel_ms = " << parallel_ms);
  INFO("timing_many_small_jobs: serial_sum   = " << serial_sum);
  INFO("timing_many_small_jobs: parallel_sum = " << parallel_sum);

  // Check correctness first
  REQUIRE(parallel_sum == Approx(serial_sum));

  if (serial_ms > 0)
  {
    double ratio = (parallel_ms > 0)
      ? double(parallel_ms) / double(serial_ms)
      : 0.0;

    INFO("timing_many_small_jobs: parallel / serial ratio = " << ratio);

    // Again: super loose bound just to catch pathological regressions.
    CHECK(ratio < 20.0);
  }
}


TEST_CASE("parallel_for: nested_serial_fallback", "[igl][parallel_for]")
{

  const int outer_loop_size = 4;
  const int inner_loop_size = 4;

  std::atomic<bool> any_inner_parallel(false);
  std::atomic<int> counter(0);

  // Outer parallel_for should use multiple threads.
  bool outer_used_parallel = igl::parallel_for(
    outer_loop_size,
    [&](int /*i*/)
    {
      bool inner_used_parallel = igl::parallel_for(
        inner_loop_size,
        [&](int /*j*/)
        {
          // Just do some work so we know the inner loop ran.
          counter.fetch_add(1, std::memory_order_relaxed);
        }
      );

      if(inner_used_parallel)
      {
        any_inner_parallel.store(true, std::memory_order_relaxed);
      }
    }
  );

  // Sanity: outer loop should be parallel when threads > 1.
  if(igl::default_num_threads() > 1)
  {
    REQUIRE(outer_used_parallel == true);
  }

  // Sanity: all iterations of both loops ran.
  REQUIRE(counter.load(std::memory_order_relaxed)
          == outer_loop_size * inner_loop_size);

  // The key assertion: inner parallel_for must fall back to serial when nested.
  REQUIRE(any_inner_parallel.load(std::memory_order_relaxed) == false);
}

#endif // IGL_PARALLEL_FOR_TIMING_TESTS
