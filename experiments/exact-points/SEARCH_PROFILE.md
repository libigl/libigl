# Isolating closest-point search costs

`profile_exact_point_search` separates tree representation/traversal from the
metric predicates used by exact nearest-triangle search. It uses the same
Armadillo input and 10,000 seeded queries as tutorial 1005. See the saved CSVs
under `results/armadillo-search-*` for timings, counters and agreement checks.

## Results

Median milliseconds for all 10,000 queries (preprocessing excluded):

| Method | One thread | Four threads |
| --- | ---: | ---: |
| Ordinary AABB | 86.03 | 22.02 |
| Ordinary Eytzinger, serialized AABB hierarchy | 99.38 | 25.25 |
| Ordinary Eytzinger, native hierarchy | 97.27 | 24.62 |
| Ordinary AABB, reconstructed Eytzinger hierarchy | 87.25 | 22.46 |
| Exact traversal, rounded comparison + rounded upper bound | 107.82 | 27.23 |
| Exact traversal, exact comparison + rounded upper bound | 132.49 | 33.76 |
| Exact traversal, rounded comparison + interval upper bound | 128.77 | 32.70 |
| Exact traversal, exact comparison + interval upper bound (copy) | 154.12 | 38.32 |
| Production exact nearest-triangle search | 150.20 | 38.14 |
| Production closest-point wrapper | 153.31 | 39.14 |
| Production cached exact-sign wrapper | 154.01 | 39.13 |
| Production uncached exact-sign wrapper | 175.70 | 44.93 |
| Uncached prototype skipping the self-comparison | 156.04 | 39.73 |
| AABB + fast winding | 97.21 | 25.06 |

**Most of the cached method's extra cost is the exact-search machinery, even
though its metric comparisons never reach expansion arithmetic on this input.**

The native ordinary Eytzinger query is 13% slower than ordinary AABB here.
With the hierarchy held fixed, Eytzinger traversal is 11–16% slower. Changing
builders while retaining the same query implementation changes time by only
about 1–2%. Thus tree representation/traversal contributes about 11 ms to the
68 ms cached-method gap over unsigned AABB, and cannot explain most of it.

The production exact search costs another 53 ms above ordinary Eytzinger.
Within the common traversal, the two comparison controls put certified leaf
ordering at an additional 25 ms and the two upper-bound controls put interval
incumbent bounds at 21 ms. These are approximate cost differences, not an
additive profiler decomposition: the benchmark copy is 2.6% slower than the
production exact search in the serial run, despite identical outputs and
traversal counts. The rounded/rounded traversal is itself about 11 ms slower
than ordinary Eytzinger; it retains validation, conservative box bounds,
different box arithmetic, and deterministic tie handling.

The four common-traversal variants all pop **3,654,186 nodes** and evaluate
**605,711 leaves**. Their incumbent updates are 72,845 (rounded comparisons)
or 72,847 (exact comparisons). All **595,711 exact comparisons** succeed in
the first `Filter` stage, as do all upper-bound feature classifications:
**zero interval fallbacks and zero expansion fallbacks**. However, the current
upper bound deliberately performs interval arithmetic on every incumbent
update, including outward rounding with `nextafter`. Ordinary native Eytzinger
evaluates 604,885 leaves, so its slightly smaller leaf count is insufficient to
explain the timing gap. Its 47,871 incumbent updates omit deterministic
lowest-index tie replacements.

Numerical closest-point reconstruction plus cached signing adds about 4 ms
to the serial exact search. The expansion audit sees 11,871 metric expansion
sums in 124 queries (1.24%) in both the closest-point and cached-sign wrappers.
Cached signing adds no metric expansions on this sample.

**Uncached signing has an additional, avoidable expansion cost.** Its vertex
edge-selection loop starts by comparing the incumbent edge with itself. The
filter cannot certify the resulting exact zero, nor can intervals, so the
comparison reaches expansions. The audit records 204,097 sums in 6,234 queries
for the production uncached wrapper. Skipping this first comparison brings it
back to the closest-point wrapper's 11,871 sums in 124 queries: it removes
192,226 sums across exactly the 6,110 closest-vertex queries. The prototype
saves about 20 ms serial / 5.2 ms with four threads, while retaining identical
signed distances and locations.

All ordinary tree variants return bit-identical squared distances on the
sample. Their triangle IDs often differ from the exact method's lowest-index
choice, but every difference is certified to be an exact distance tie. The
exact traversal copy matches every production triangle ID. Cached, uncached
and skip-self signed results match bit for bit. This workload uses uniform
padded-box samples; near-surface or adversarial inputs can have different
feature distributions and fallback rates.

The initial measurements motivated follow-up controls for repeated feature
classification and pruning-bound recomputation on ties, described below.
A cheap certified floating-point upper-bound filter before the interval path
is another possible optimization that has not been measured here.

## Shared-feature ties

The production comparator already short-circuits exact distance ties when both
closest-feature classifications identify the same vertex coordinates or the
same edge endpoints, irrespective of edge orientation. It performs this check
before constructing distance fractions. The expensive expansion fallback seen
in uncached vertex signing does **not** occur in these traversal comparisons.

In the Armadillo traversal, 67,779 of 595,711 comparisons (11.38%) are ties:
63,629 at shared vertices and 4,150 at shared open edges. Every one is recognized
by the existing shared-feature check. But the two feature classifications are
repeated for every comparison, and 25,117 ties select a lower triangle ID and
recompute the incumbent's interval upper bound.

The follow-up prototypes separate three changes:

- Retain the previous upper bound on an exactly equal distance, including when
  the deterministic winning triangle ID changes. This is safe because a bound
  on the old distance is a bound on the equal new distance. It reduces interval
  bound evaluations from 72,847 to 47,730 without changing the visited nodes.
- Save only the incumbent's feature dimension and local feature ID during the
  query. Classify each new candidate once, reuse these classifications for the
  comparison and upper bound, and retain the full filtered/exact arithmetic
  chain for comparing different features. No expansions survive a comparison.
- Recognize shared features either by coordinate equality (as in production)
  or by mesh vertex indices. Both versions cache the same feature metadata.
  Comparing them distinguishes the value of topology from avoiding repeated
  classification. With distinct indices at coincident coordinates, the topology
  version can miss the shortcut but still resolves the tie through exact math.

All versions preserve the lowest triangle index on a tie. The instrumented
prototypes visit the same 605,711 leaves and return the same 10,000 triangle IDs
as production. An additional 180 checks compare production and the five search
prototypes against brute-force exact scans, using on-surface points, near ties,
reversed face orderings and coincident vertices with different indices.
The saved follow-up runs are `results/armadillo-ties-*`.

| Search variant | One-thread median, ms |
| --- | ---: |
| Production exact search | 152.07 |
| Benchmark copy, original policy | 151.23 |
| Retain upper bound on ties | 145.01 |
| Cache features, compare coordinates | 133.78 |
| Cache features, compare vertex indices | 134.83 |
| Cache features + retain tie bound, coordinates | 127.52 |
| Cache features + retain tie bound, vertex indices | 126.99 |

The combined prototype is about 16% faster than production exact search on this
sample. Index equality has no clear speed advantage over coordinate equality
once both versions cache the classifications. The savings come from reusing
information: only 605,711 candidate feature classifications are needed instead
of the original 1,264,269 classifications (two per comparison plus one per
incumbent update, with no classification fallbacks on this sample).

Tie-bound reuse alone saves about 6 ms in this run; an earlier run measured
8 ms. Feature caching also benefits non-tied comparisons, so its entire saving
must not be attributed to ties. All search prototypes use zero metric expansion
sums on Armadillo, independently confirmed by the audit executable.

The follow-up four-thread run passed all accuracy checks but overlapped another
CPU-heavy process sharing this container's four-core quota. Its timing ranges
are wide and several medians are distorted; it is preserved as
`armadillo-ties-4-threads-contended.csv` and is **excluded from performance
conclusions**. The initial one-/four-thread tree-versus-predicate controls above
were collected earlier. These optimizations remain experiment-local prototypes;
the production library has not been modified in this investigation.

## Controls

All methods prepare their data before timing, reuse output allocations, and
run scalar queries under the same `igl::parallel_for` with threshold 10,000,
explicit thread count, and worker exception handling. There are no nested
query pools. Eigen uses one thread. Each method gets one warmup batch, followed
by nine measured batches with rotating method order; the table reports medians.
Loading, preprocessing, correctness checks and counters are excluded.

The ordinary AABB and Eytzinger methods both use
`point_simplex_squared_distance<3>` at leaves and return squared distance,
triangle index and closest point. We test both native builders, serialize the
AABB hierarchy for Eytzinger traversal, and reconstruct an AABB from the native
Eytzinger hierarchy. Those conversions preserve the partition and boxes. They
control for builder differences, but their timings still combine storage
layout and traversal algorithm: AABB recurses and checks child containment;
Eytzinger uses an explicit active-node vector and different box arithmetic.

`search_ablation.h` copies the production exact traversal and independently
switches its leaf comparison and incumbent upper bound. All four variants keep
its box lower bound, validation, explicit stack, traversal order and tie rule.
The exact/exact copy must return exactly the production triangle ID for every
query. Rounded comparison uses the ordinary triangle-distance kernel. Rounded
upper bounds use that kernel only when the incumbent changes (or reuse its
result if the comparison was also rounded).

**The rounded ablations are cost controls, not robust query implementations.**
Their pruning bounds and/or candidate ordering are inexact. Sample agreement
does not establish their correctness on other inputs. The production exact
search, closest-point wrapper, cached/uncached signed wrappers and fast-winding
query are also timed directly. The fast-winding method preserves libigl's
continuous `sqrt(squared_distance) * (1 - 2*abs(w))` formula.

The untimed traversal audit counts popped nodes, leaf evaluations, incumbent
updates and Filter → Interval → Expansion fallbacks. A separate executable,
`audit_exact_point_search`, interposes `predicates_expansion_sum` to count metric
expansion operations and queries using them. It runs serially; the timed
executable has no arithmetic interposition. This counter does not count
Shewchuk's internal orient3d arithmetic.

`prototype_uncached_skip_self` copies only the local vertex sign logic, starting
its edge loop after the incumbent edge. It retains the exact predicates for all
remaining comparisons. Its full signed-distance outputs must equal both
production cached and uncached outputs bit for bit on the sample. This is an
experiment-local prototype; the production signing implementation is unchanged.

## Reproduce

Configure `experiments/exact-points` using the existing dependency directories
as described in [README.md](README.md). The input is the retained tutorial
`armadillo.obj`; no additional dataset is needed.

```sh
cmake --build build-exact-points \
  --target profile_exact_point_search audit_exact_point_search -j 4
taskset -c 24 build-exact-points/profile_exact_point_search \
  tutorial/data/armadillo.obj 10000 9 1
taskset -c 24-27 build-exact-points/profile_exact_point_search \
  tutorial/data/armadillo.obj 10000 9 4
taskset -c 24 build-exact-points/audit_exact_point_search \
  tutorial/data/armadillo.obj 10000 1 1
```

Use affinity masks available on the machine. Run each thread count in a new
process because libigl caches its thread count. The expansion audit requires
Linux and a GNU-compatible linker; timing controls do not require interposition.
Recorded runs use GCC 13.3, `-O3 -DNDEBUG -fno-fast-math -ffp-contract=off`, an
AMD EPYC 9B45, and no concurrent compilation. Input revision, hash and sampling
metadata are in `results/validation-armadillo-profile.json`.
