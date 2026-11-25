# Certified closest-point filters

This round adds three optimizations to `7264646f`:

1. Obtain an upper bound from the filtered distance numerator and denominator,
   using three outward-rounding operations; retain the interval fallback.
2. Keep the incumbent's ordinary filtered fraction and absolute permanents,
   avoiding their recomputation at every leaf comparison.
3. Reuse dot products between closest-feature classification and distance
   evaluation. Uncertain classifications still use the existing predicate.

All expansions remain local to fallback calls. The public API, smallest-triangle
id rule, coordinate-based shared-feature shortcuts and sign calculations are
unchanged. This work follows the feature caching, exact-tie bound reuse and
vertex self-comparison fixes documented in [INTEGRATION.md](INTEGRATION.md).

## Isolated search costs

Armadillo, 43,243 vertices and 86,482 faces; 10,000 identical padded-box queries,
seed 896; GCC 13.3, Release `-O3 -DNDEBUG`, strict FP. One warm-up then nine
rotated repetitions. One thread pinned to CPU 24; four threads pinned to 24–27.
The same `parallel_for` policy and threshold of 10,000 apply to every variant.
No builds or other observed CPU-heavy processes ran alongside the measurements.

Median milliseconds per complete batch (excluding setup and point construction):

| Variant | 1 thread | 4 threads |
| --- | ---: | ---: |
| Previous traversal copy | 127.20 | 33.54 |
| Filtered upper bounds only | 113.26 | 29.84 |
| Reused distance only | 124.74 | 32.87 |
| Reused distance + filtered upper | 111.09 | 29.56 |
| Fused classification + reuse + filtered upper | **107.90** | **28.76** |
| Full distance enclosures + reuse | 123.18 | 32.60 |
| Full enclosures + fused classification + reuse | 119.80 | 32.48 |
| Previous production implementation | 126.19 | 33.71 |

The retained combination saves 15.2% against the serial baseline copy. The
smaller changes are incremental improvements; their percentages do not add.
The independent production control also checks that copying the traversal into
the harness has not substantially changed its cost. The preliminary unpinned
run gave the same ranking and is retained separately.

Computing a full distance enclosure for every leaf eliminates all non-topological
cross-product comparisons on this sample, but its endpoint rounding and division
costs exceed those comparisons. It remains a benchmark control, not a production
path. Merely reducing the number of predicate calls is insufficient: their
common filtered path is already cheap.

Every variant visited exactly 3,654,186 nodes and 605,711 leaves, made 47,730
strict incumbent updates, and returned identical triangle ids. The winning
variant certified all feature classifications and all upper bounds in its first
filter; no interval or expansion distance-comparison fallbacks occurred. Its
527,932 comparisons excluding shared-feature ties use the cached cross-product
filter. The former traversal performed 47,730 interval upper bounds. Counts come
from separate untimed instrumented copies, with counter branches compiled out
of the timed variants.

Raw results: `results/armadillo-filters-ablation-{serial,four,unpinned}.csv`.

## Why the bounds remain conservative

For each polynomial, `Filter` stores its rounded value and a positive absolute
permanent. The existing `256*epsilon*permanent` error bound applies to both the
old and fused formulas. Using the recurrence in [README.md](README.md), the five
base dot products have rounding weight 5. Subtracted edge projections have
weight at most 7; the reused second-edge length has weight 7; its numerator has
weight 16. Face numerator/denominator weights remain 17/11, and the largest
comparison has weight 30. The existing bound therefore still covers all of them.

For certified nonnegative N and positive D, the upper bound is
`up(up(N + error_N) / down(D - error_D))`. Each endpoint and the quotient round
outward. An exactly zero permanent and numerator permit upper bound zero.
Nonfinite or underflowed error estimates, uncertified numerator/denominator,
and nonfinite quotients trigger the original interval upper-bound calculation.
This preserves the documented finite-range and strict-FP contract; it does not
extend the predicates to arbitrary exponents.

The fused classifier uses exactly the existing vertex-first, edge-second region
tests, algebraically reusing their dot products. An uncertain sign immediately
returns to the existing classifier. Refinement and exact comparison use original
coordinates, never a rounded projection or cached approximate feature decision.

After integration, a separate run compared production with the retained controls
in the same executable: **125.38 ms baseline copy, 108.16 ms production**
(13.7% reduction), with 106.31 ms for the winning template prototype. All outputs
still match. This small prototype/production difference illustrates code-generation
sensitivity; the production measurement is the relevant implementation result.
See `results/armadillo-filters-production-serial.csv`.

## Complete signed-distance comparison

The before binary is from `7264646f`; the after binary differs in the two private
point-kernel files. Tutorial source, compiler options, dependency and mesh are
identical. For each version and thread count, three processes each warmed all
methods once and measured nine rotated repetitions. Order of before/after
processes was reversed between pairs. Tables show the median of the three
per-process medians. All runs, including noisier ones, are retained.

Query milliseconds for 10,000 points with preprocessing already done:

| Method | 1 thread before | 1 thread after | 4 threads before | 4 threads after |
| --- | ---: | ---: | ---: | ---: |
| Unsigned squared distance | 109.80 | 112.32 | 26.52 | 27.31 |
| Pseudonormal signed distance | 99.41 | 101.21 | 24.23 | 24.97 |
| Fast winding signed distance | 117.03 | 119.73 | 28.37 | 29.50 |
| Exact sign, uncached | 139.38 | 131.70 | 39.21 | 32.34 |
| Exact sign, cached | 145.21 | 123.99 | 37.31 | 30.99 |

Four-thread exact-sign time decreased 17.5% uncached and 16.9% cached. Serial
reductions are 5.5% and 14.6%, respectively. These complete-wrapper runs varied
more than the isolated search runs: serial per-process cached baseline medians
were 130.48, 145.21 and 145.81 ms, while the unchanged unsigned control ranged
89.58–109.93 ms. Methods' independent medians can even put cached baseline time
above uncached time. The source-unchanged controls shift 2–4% between versions;
do not interpret the table as a precise universal speedup or a precise remaining
gap against fast winding. The isolated controls provide the stronger evidence
for which changes help.

Median full-call milliseconds including preprocessing:

| Method | 1 thread before | 1 thread after | 4 threads before | 4 threads after |
| --- | ---: | ---: | ---: | ---: |
| Unsigned squared distance | 174.45 | 156.51 | 78.73 | 82.93 |
| Pseudonormal signed distance | 223.09 | 222.19 | 119.51 | 120.55 |
| Fast winding signed distance | 209.80 | 222.97 | 97.46 | 100.47 |
| Exact sign, uncached | 235.13 | 198.03 | 109.60 | 105.93 |
| Exact sign, cached | 392.14 | 408.73 | 279.22 | 274.83 |

Setup algorithms are unchanged. Full-call results have wider variation and do
not demonstrate an improvement for the cached serial call: its reported median
increases 4.2%, while per-process baseline medians span 387.46–474.33 ms and
optimized medians span 401.05–434.87 ms. Setup dominates this path; these runs
do not isolate a setup regression. Exposure caching still saves only the local
sign stage, and its preprocessing is expensive for a single 10,000-point batch.
The query improvement comes from cheaper certified arithmetic, not fewer nodes
or fewer exact expansion fallbacks.

Raw complete-wrapper results are
`results/armadillo-filters-{baseline,optimized}-{serial,four}*.csv`.
Every run preserves identical cached/uncached exact outputs and prepared/full-call
results. All sampled signs agree with the pseudonormal and winding methods.

## Validation and reproduction

- Release: 88,003 assertions in 24 test cases; all three CTest targets pass,
  including explicit static-library linkage and the allocation/path audit.
  The same tests pass with a 1 MiB stack.
- New tests check fused classifications against the independent bounded-integer
  projection oracle and 3,000 expansion-reference cases, including broad
  exponents, exact coincidences, degeneracy and one-ulp perturbations. They
  certify the numerator, denominator, upper bound and reused comparison result.
  Dedicated cases require cancellation/zero refinement and rejection of
  nonfinite or underflowed error estimates.
- The production path regression fails against `7264646f` (45 rounding calls)
  and passes now (3), for both a single leaf and an exact-tie incumbent replacement.
  Leaf predicates still allocate nothing; the cancellation cases still perform
  18 expansion sums for ordering and 35 for numerical construction.
- ASan/UBSan with leak detection: all three CTest targets pass. The full
  10,000-point, four-thread Armadillo tutorial also passes, including all output
  agreement checks. Sanitizer timings are not performance measurements.
- The standalone final production search agrees on every triangle id with the
  prior traversal copy. The interval-only ablation patch passes `git apply --check`
  and disables the new filtered entry points as well as the original ones.

The latest changes touch only the point kernel. The ray and dependency validation
from [INTEGRATION.md](INTEGRATION.md) remains applicable; those paths are unchanged.
Machine-readable configuration and results are in `results/validation-filters.json`.

Use the existing-dependency configure command in [README.md](README.md), including
`TUTORIAL_DATA_DIR` for the complete example, then:

```sh
cmake --build build-exact-points --target profile_point_filters test_exact_points static_smoke audit_allocations 1008_SignedDistanceProfile -j 4
ctest --test-dir build-exact-points --output-on-failure
(ulimit -s 1024; IGL_NUM_THREADS=4 build-exact-points/test_exact_points)
taskset -c 24 build-exact-points/profile_point_filters tutorial/data/armadillo.obj 10000 9 1
taskset -c 24-27 build-exact-points/profile_point_filters tutorial/data/armadillo.obj 10000 9 4
taskset -c 24 build-exact-points/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 9 1
taskset -c 24-27 build-exact-points/1008_SignedDistanceProfile tutorial/data/armadillo.obj 10000 9 4
```

For paired complete-wrapper timings, build one binary at `7264646f` and one with
these changes using identical compiler/dependency configuration. Preserve the
first executable before rebuilding. The isolated harness includes the prior
traversal policy and all controls directly, with no need for a second checkout.
Comparison-stage counters apply to the variants caching fractions; fused-feature
fallback counters apply to fused variants. Counter zero in other variants means
that stage was not instrumented, not that no work was done.
