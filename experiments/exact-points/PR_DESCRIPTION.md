# Use Shewchuk expansions for exact ray and closest-point mesh queries

This extends #2549 with exact closest-feature selection, closest-triangle
ordering and local inside/outside/on-surface tests. The signed-distance wrapper
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
- Cached feature classifications during traversal, upper-bound reuse on exact
  ties, and removal of the incumbent edge's self-comparison in vertex signing.
  Different features still use the complete filtered/exact comparison path.
- Headless tutorial `1008_SignedDistanceProfile`, comparing unsigned,
  pseudonormal, fast-winding and cached/uncached exact-sign distances with the
  same query parallelism. The branch is rebased onto current main `cb3d2940`.

## Performance

Armadillo, 86,482 triangles and 10,000 fixed padded-box queries, GCC 13.3 Release,
strict FP, EPYC 9B45. Times exclude preprocessing. Before/after binaries share
the same rebased source except for the three optimized point-kernel files.

| Method | 1 thread before → after | 4 threads before → after |
| --- | ---: | ---: |
| Exact sign, uncached | 183.81 → 137.02 ms | 51.02 → 37.11 ms |
| Exact sign, cached | 154.65 → 130.39 ms | 42.99 → 35.07 ms |
| Fast winding, control | 94.40 → 94.28 ms | 26.47 → 26.03 ms |

Serial query time drops 25.5% uncached and 15.7% cached. Cached exact-sign queries
remain about 38% slower than fast winding serially on this sample. Exposure
preprocessing remains expensive for a single batch; caching now saves about
5% of query time. Four-thread values aggregate three runs and show more variance.
[Methodology, full-call times, raw results and work counters](INTEGRATION.md)
are included; [API contracts and usage](README.md) describe the supported inputs.

## Validation

- Point tests: 56,844 assertions in 22 cases; all three Release and ASan/UBSan
  CTest targets pass, including static linkage and allocation/path audits.
  The Release suite also passes with a 1 MiB stack.
- Both targeted optimization-path regressions fail before the changes and pass
  afterward. Cancellation tests still exercise expansion fallbacks.
- All 10,000 Armadillo exact triangle ids agree with the original-policy search;
  cached/uncached signed outputs are bit identical. The full four-thread example
  passes ASan/UBSan with leak detection. The standard tutorial target also builds.
- Rays: 22,295 assertions in 13 cases, all five CTest targets, and comparison to
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

Prepared by Codex.
