# Use Shewchuk expansions for exact ray and closest-point mesh queries

Supersedes the earlier draft #2549, carrying forward its ray implementation and
adding exact closest-feature selection, closest-triangle ordering and local
inside/outside/on-surface tests. The signed-distance wrapper
returns an approximate numerical distance together with exact triangle/feature
ids and an exact location label. It supports both single points and batches
against a caller-owned Eytzinger tree and mesh data.

The branch includes the ray implementation from #2549 and moves its exact
operations into `igl::predicates`. Public ray-triangle classification and hit
ordering functions use Shewchuk's existing arithmetic with local stack buffers;
ray-box traversal remains in core. A small CMake dependency patch exposes the
existing expansion primitives and prevents out-of-bounds lookahead reads.

For closest-point queries, the branch adds:

- Exact vertex/open-edge/open-face classification and rational squared-distance
  comparison, including degenerate triangles and deterministic smallest-id ties.
- Conservative Eytzinger traversal, local mesh signing using array triangle
  adjacency, and optional query-independent edge/vertex exposure caching.
- `point_mesh_closest_point` and `approximate_signed_distance_exact_sign`,
  retaining approximate point coordinates, barycentrics and distance magnitudes
  separately from the exact selection and sign information.
- Cached feature classifications and filtered distance fractions during traversal,
  fused feature/distance evaluation, cheap certified upper bounds, upper-bound
  reuse on exact ties, and no incumbent self-comparison in vertex signing.
  Uncertain decisions retain the interval/expansion refinement path.
- Headless tutorial `1008_SignedDistanceProfile`, comparing unsigned,
  pseudonormal, fast-winding and cached/uncached exact-sign distances with the
  same query parallelism.

## Performance

Armadillo, 86,482 triangles and 10,000 identical padded-box query points,
GCC 13.3 Release, strict floating-point operations, EPYC 9B45. The following
complete signed-distance query times are directly comparable: **four threads,
preprocessing excluded, all distance and signing work included**, and the same
parallel policy for every method.

| Method | Median query time for 10,000 points |
| --- | ---: |
| Fast winding signed distance | 29.50 ms |
| Approximate distance with exact sign, cached | 30.99 ms |
| Approximate distance with exact sign, uncached | 32.34 ms |

Each value is the median of three process medians, with one warm-up and nine
rotated repetitions per process. Cached exact signing is about 5% slower than
fast winding in these runs. Full-wrapper timings vary, so this is a workload
measurement rather than a precise general speedup claim. Exposure-cache setup
is excluded here; setup dominates cached full-call time and those serial
measurements showed no improvement from the latest filters.

Separate **one-thread closest-triangle search** controls measure 108.16 ms for
the integrated implementation versus 125.38 ms for the preceding traversal
policy (13.7% less). Cheaper certified upper bounds save most of this cost;
reusing filtered fractions and fusing classification add smaller improvements.
All variants visit the same leaves and select identical triangle ids.

These measurements use the tested base `cb3d2940`. Current main has since
advanced, including a fast-winding second-moment correction; that correction
is not included in these benchmark numbers. [All paired timings, methodology,
bounds proof and counters](FILTERS.md) and [API contracts and usage](README.md)
are included.

## Validation

- Point tests: 88,003 assertions in 24 cases; all three Release and ASan/UBSan
  CTest targets pass, including static linkage and allocation/path audits.
  The Release suite also passes with a 1 MiB stack.
- The new upper-bound path regression fails against the preceding implementation
  (45 rounding calls) and passes after the change (3). Cancellation still reaches
  expansion fallbacks. New certificates are checked against full expansions,
  including broad exponents, one-ulp perturbations and degenerate triangles.
- All 10,000 Armadillo exact triangle ids agree with the original-policy search;
  cached/uncached signed outputs are bit identical. The full four-thread example
  passes ASan/UBSan with leak detection. The standard tutorial target also builds.
- Unchanged ray paths previously passed 22,295 assertions in 13 cases, all five
  CTest targets, and comparison to
  #2549 on 125,000 classifications, 75,031 hit comparisons and 1,200 mesh queries
  with zero mismatches in the tested domain.
- The dependency patch applies cleanly to the pinned upstream source and is
  idempotent. The tested patched source matches the companion dependency.

Exactness assumes the documented arithmetic range and strict floating-point
operations. Signing requires an embedded, closed, consistently outward-oriented
triangle manifold with nondegenerate faces; cache preprocessing does not certify
self-intersection freedom or global shell orientation. Numerical distances and
constructed points are approximate. Existing ray-box culling and approximate-hit
limitations inherited from #2549 are documented in
[the ray report](../exact-rays/README.md#accuracy-and-limitations).

#### Checklist

- [ ] All changes meet [libigl style-guidelines](https://libigl.github.io/style-guidelines/).
- [x] Adds new .cpp file.
- [x] Adds corresponding unit test.
- [ ] This is a minor change.

Prepared by Codex (AI agent).
