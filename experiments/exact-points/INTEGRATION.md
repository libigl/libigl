# Production closest-point optimizations

This report records the first integration round. See [FILTERS.md](FILTERS.md)
for subsequent certified-bound and distance-reuse optimizations and current timings.

The closest-point traversal now classifies each candidate once, keeps the
incumbent feature dimension and local id, and reuses its conservative upper
bound when a smaller triangle id replaces it on an exact distance tie.
Shared-vertex and shared-edge shortcuts still compare coordinates, including
coincident vertices with distinct indices and reversed edge orientations.
Other distance comparisons retain the Filter → Interval → Expansion path.

Uncached vertex signing starts after the incumbent edge in the validated
one-ring, avoiding an exact-zero self-comparison. Other edge comparisons still
use the exact metric predicate. No expansion or constructed exact point is
retained across calls, and the public interfaces and smallest-face-id tie rule
are unchanged.

The branch is rebased onto main `cb3d2940`, including the upstream thread pool
and pseudonormal changes. The headless example is now tutorial
`1008_SignedDistanceProfile`; main already uses 1005–1007. The before/after
binaries use this same base and example. Only the three production point-kernel
files differ: the before version is from `7344d884`, the after version from
`ad7d7fe9`. This isolates these optimizations from the rebase.

## Armadillo timings

43,243 vertices, 86,482 triangles, the same 10,000 uniform padded-box queries
(seed 896), GCC 13.3, Release `-O3 -DNDEBUG -fno-fast-math -ffp-contract=off`,
AMD EPYC 9B45, a four-core container quota. All five prepared methods use the
same `igl::parallel_for` threshold (10,000), worker policy, exception handling,
output reuse and thread count; Eigen uses one thread. Setup is excluded from
query times. Each process warms every method once, then rotates method order
through nine measured repetitions. No compilation ran during timings.

Serial runs were pinned to CPU 24. Four-thread runs were pinned to CPUs 24–27
and repeated in three processes per version because the first run had wider
ranges. All runs are retained; the four-thread columns are the median of the
three per-process medians, with before/after order reversed between runs.
There was no observed competing CPU-heavy process, but four-thread timings and
especially full-call setup times still vary. Treat their differences as workload
measurements, not precise general speedup guarantees.

Median query time for the complete 10,000-point batch:

| Method | 1 thread before | 1 thread after | 4 threads before | 4 threads after |
| --- | ---: | ---: | ---: | ---: |
| Unsigned squared distance | 88.51 ms | 88.46 ms | 25.16 ms | 24.21 ms |
| Pseudonormal signed distance | 83.65 ms | 82.25 ms | 24.03 ms | 22.02 ms |
| Fast winding signed distance | 94.40 ms | 94.28 ms | 26.47 ms | 26.03 ms |
| Exact sign, uncached | 183.81 ms | 137.02 ms | 51.02 ms | 37.11 ms |
| Exact sign, cached | 154.65 ms | 130.39 ms | 42.99 ms | 35.07 ms |

Serial exact-sign query time decreased by **25.5% uncached** and **15.7% cached**.
The four-thread median reductions are **27.3%** and **18.4%**, respectively.
Cached exact-sign queries remain about **38% slower than fast winding** serially
and **35% slower with four threads** on this sample. Exposure caching now saves
about 5% of query time; its roughly 170 ms serial preprocessing is worthwhile
only with substantial reuse for this distribution.

Median full-call time, including tree/adjacency/cache or normal/winding setup:

| Method | 1 thread before | 1 thread after | 4 threads before | 4 threads after |
| --- | ---: | ---: | ---: | ---: |
| Unsigned squared distance | 139.95 ms | 139.93 ms | 79.75 ms | 80.77 ms |
| Pseudonormal signed distance | 179.39 ms | 178.15 ms | 119.03 ms | 122.38 ms |
| Fast winding signed distance | 186.97 ms | 187.30 ms | 96.39 ms | 98.49 ms |
| Exact sign, uncached | 257.62 ms | 212.77 ms | 123.88 ms | 116.58 ms |
| Exact sign, cached | 401.07 ms | 381.98 ms | 290.72 ms | 289.57 ms |

Setup is unaffected by these optimizations. In particular, cached full-call
cost is dominated by serial exposure preprocessing; the four-thread total-time
difference is within the observed variability. The raw runs include all
minima, maxima, component setup timings, checksums and accuracy reports under
`results/armadillo-integrated-{baseline,optimized}-*`. The one-thread result has
one process per version; the four-thread result includes the unnumbered file
and `-run2` / `-run3` files. These are final runs on the rebased branch.

## Accuracy and removed work

The Armadillo audit preserves all 10,000 exact triangle ids relative to the
original-policy traversal copy and the optimized prototypes. It also passes
180 brute-force comparisons for shared features, coincident coordinates with
different vertex indices, reversed faces, exact surface points, near ties and
large query coordinates. Cached and uncached signed outputs are bit identical.
All sampled signs agree with the other signed methods; prepared and convenience
API distances agree bit for bit. These samples do not replace the arithmetic
contract or establish correctness for arbitrary exponent ranges.

The metric expansion audit now reports:

| Path | Expansion sums | Queries using expansions |
| --- | ---: | ---: |
| Exact closest-triangle search | 0 | 0 |
| Numerical closest-point construction | 11,871 | 124 |
| Cached signed wrapper | 11,871 | 124 |
| Uncached signed wrapper | 11,871 | 124 |

The previous uncached path performed 204,097 sums in 6,234 queries. Removing
self-comparisons eliminates 192,226 sums across the 6,110 closest-vertex queries.
These counters cover metric expansion wrappers, not orient3d's internal work.
The matching traversal prototypes need 605,711 feature classifications instead
of 1,264,269, and 47,730 interval upper bounds instead of 72,847, visiting the
same 3,654,186 nodes and 605,711 leaves. Production query output matches them;
the node/classification counters belong to the instrumented copies.

Two focused production-path tests were run before and after the fixes. For a
duplicate-triangle exact tie, a single query and a two-leaf query now both make
45 outward-rounding calls; the old two-leaf query made 90. A tetrahedron vertex
sign with distinct edge projections now makes zero expansion sums instead of
16. The unchanged cancellation cases still reach expansions (18 sums for leaf
ordering and 35 for numerical construction), and predicates allocate nothing.

## Validation

- Release: all three point CTest targets passed, with 56,844 assertions in
  22 cases, including 1,000 comparisons using cached classifications against
  full expansions and 30 tree regressions for ties and nearly coincident faces.
  The same suite passed with a 1 MiB stack.
- AddressSanitizer and UndefinedBehaviorSanitizer, with leak detection: all
  three point CTest targets and the complete 10,000-point, four-thread Armadillo
  example passed. Static-library APIs and worker exception propagation are covered.
- Rays: all five CTest targets passed; the complete unit suite including the
  cube test passed 22,295 assertions in 13 cases. Comparison to the original
  PR checked 125,000 classifications, 75,031 hit comparisons and 1,200 mesh
  queries with zero mismatches in the supported test domain.
- The normal headless tutorial CMake target builds and passes a 16-point smoke
  run in addition to the full experiment-harness runs.
- The dependency patch applies to the pinned upstream source, matches the
  tested companion dependency byte for byte, and tolerates repeated application.

Machine-readable metadata is in `results/validation-integration.json`. Existing
[ray limitations](../exact-rays/README.md#accuracy-and-limitations), including
core ray-box failures and approximate hit output, are outside these point-query
optimizations. Exact metric results retain their documented finite-range and
strict floating-point assumptions; signed distances themselves are approximate.

## Reproduce

Configure and build with the existing dependencies as in [README.md](README.md).

```sh
cmake --build build-exact-points -j 4
ctest --test-dir build-exact-points --output-on-failure
(ulimit -s 1024; IGL_NUM_THREADS=4 build-exact-points/test_exact_points)
taskset -c 24 build-exact-points/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 9 1
taskset -c 24-27 build-exact-points/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 9 4
taskset -c 24 build-exact-points/audit_exact_point_search tutorial/data/armadillo.obj 10000 1 1
```

The saved before binary was built using the three point-kernel files from the
parent of the optimization commit, with all other sources and compiler options
held fixed. Restore the optimized files and rebuild before running the path
audit: those tests intentionally fail with the old code. The interval-only
ablation patch was updated and passes `git apply --check`; historical
interval-only timings predate these traversal changes.
