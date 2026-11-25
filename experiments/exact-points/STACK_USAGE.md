# Exact metric predicates on small worker stacks

The baseline is `caafc200`. Its static closest-point construction exceeded a
512 KiB macOS worker stack. This also reproduces with GCC 13.3 on Linux:
the new `exact_point: 512_KiB_worker_stack` test exits with SIGSEGV (139) against
the baseline static library. The existing allocation audit likewise crashes
with `ulimit -s 512`. The earlier 1 MiB test did not cover this limit.

## Change

Expansion storage now follows expression bounds. An input has one component,
a difference of input coordinates has at most two, a dot product of differences
has at most 24, and a cross-product component has at most 16. Addition and
subtraction propagate `N1+N2`; multiplication propagates `2*N1*N2`. Storage is
capped at the existing 2098-component binary64 bound. Temporary sum/scale
buffers still hold the full possible output before its length is checked.
Conversions can only widen capacity, and arithmetic order is unchanged.

The exact construction is called through a non-inlined helper, keeping its
large return slot off the filtered query's frame. All arithmetic storage is
still local to the predicate; no thread-local storage or heap allocation was
introduced, and callers do not need to enlarge their thread stacks.

GCC 13.3 `-O3 -DNDEBUG -fno-fast-math -ffp-contract=off -fstack-usage` reports
these maximum individual frames across the static point API translation units:

| Function | Baseline bytes | Fixed bytes |
| --- | ---: | ---: |
| `point_mesh_closest_point` | 436,928 | 304 |
| `construct<Expansion>` | 352,864 | 118,432 |
| `construct_exact` (new cold entry) | — | 134,384 |
| `feature<Expansion>` | 352,848 | 64,208 |
| `distance<Expansion>` | 386,384 | 61,232 |
| `compare<Expansion>` | 201,584 | 201,584 |
| `perpendicular_compare<Expansion>` | 386,384 | 123,712 |
| `exact_mesh_feature_exposure` | 336,592 | 2,112 |

The frame measurements are saved in `results/stack-frames.json`. These are
compiler-specific frames, not whole-call stack peaks. The worker
test verifies nested execution with an explicitly reserved 512 KiB stack,
using pthread attributes on POSIX and `CreateThread` on Windows. It covers
distance refinement, exact ties, feature classification, exposure preprocessing,
near-plane construction, on-surface edges, and local/cached signing. The
allocation audit independently confirms that the near-plane construction
executes 35 expansion sums while its feature classification needs none.

## Validation

- Release: 88,007 assertions in 25 cases; all five experiment CTest targets pass.
- ASan/UBSan with leak detection: all five CTest targets pass, including the
  512 KiB worker test in both header-only and static builds.
- The full 10,000-point, four-thread Armadillo tutorial also passes ASan/UBSan
  with leak detection and zero sign mismatches.
- The allocation audit passes with a 512 KiB process stack. Leaf arithmetic
  still allocates nothing; each wrapper query retains its one traversal-vector
  allocation, with none added by exact construction.
- All 12 Armadillo benchmark runs have identical query/full-call checksums,
  agreement reports, and feature/location counts. Each run checks all 10,000
  cached/uncached results and reports zero sign mismatches.

These measurements were made on Linux x86-64 with GCC. Native macOS and Windows
coverage is provided by the new regression in the regular predicate test suite;
those CI results must be checked separately.

## Performance

Same tutorial, mesh and seed as the previous reports: 86,482 Armadillo faces,
10,000 points, seven rotated-order timing repetitions per process, three
processes per variant/thread count. Baseline and fixed binaries were alternated
between processes. One thread was pinned to CPU 24; four to CPUs 24–27.
Eigen used one thread and every method used the same libigl parallel policy.
No builds ran during measurement. Values below are medians of the three
process medians, in milliseconds.

| Measurement | Threads | Baseline | Fixed |
| --- | ---: | ---: | ---: |
| Fast winding, prepared query | 1 | 96.19 | 94.59 |
| Exact sign, uncached prepared query | 1 | 118.59 | 116.40 |
| Exact sign, cached prepared query | 1 | 112.03 | 109.69 |
| Fast winding, prepared query | 4 | 24.66 | 24.27 |
| Exact sign, uncached prepared query | 4 | 30.22 | 29.69 |
| Exact sign, cached prepared query | 4 | 28.66 | 27.88 |
| Exact sign, cached setup + query | 1 | 364.42 | 307.16 |
| Exact sign, cached setup + query | 4 | 269.90 | 214.83 |

There is no observed query-time regression. The unchanged fast-winding control
also improves by about 1.6–1.7%, so the small query gains should not be treated
as a reliable speedup. The cached full-call improvement (16–20%) is larger:
exposure preprocessing also benefits from the smaller exact arithmetic frames.
The synthetic mesh microbenchmark reduces closest-point-wrapper time from
1.80 to 1.50 us for 320 faces and 5.93 to 5.44 us for 5,120 faces. Its exposure
preprocessing is about 3x faster; ordinary leaf filters are essentially unchanged.

Raw runs are in `results/armadillo-stack-{baseline,fixed}-{1,4}-thread-run{1,2,3}.csv`;
`results/stack-timings-summary.json` retains every process median. The synthetic
measurements are in `results/timings-stack-{baseline,fixed}.csv`.

Build with the dependency overrides in [README.md](README.md#reproduce), adding
`-DCMAKE_CXX_FLAGS=-fstack-usage` to inspect GCC frames. Run `static_stack` and
`test_exact_points 'predicates::exact_point: 512_KiB_worker_stack'` for the explicit
worker test. For Armadillo timings, build `1008_SignedDistanceProfile` with
`TUTORIAL_DATA_DIR` set and invoke it as:

```sh
taskset -c 24 build/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 7 1
taskset -c 24-27 build/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 7 4
```
