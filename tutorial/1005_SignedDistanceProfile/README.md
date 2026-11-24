# Signed-distance profiling

A headless example comparing five methods on the tutorial `armadillo.obj`
(43,243 vertices, 86,482 triangles):

- `point_mesh_squared_distance` (unsigned, returning squared distances).
- `signed_distance` with angle-weighted pseudonormals.
- `signed_distance` with fast winding number (order 2, accuracy scale 2).
- `predicates::approximate_signed_distance_exact_sign` using triangle adjacency
  to compute local signs without an exposure cache.
- The same function with a precomputed `MeshFeatureExposure` cache.

By default it samples 10,000 points uniformly in the mesh's axis-aligned bounding
box, padded by 10% of each axis extent on each side. The seed is fixed at 896.
Every method receives the identical double-precision points and mesh. The
example checks manifoldness, boundary edges, consistent adjacency, degenerate
faces, positive signed volume and self-intersections before profiling queries.
Custom meshes must be outward-oriented, including the appropriate orientation
of any cavity shells; a positive total volume alone cannot verify shell nesting.

## Build and run

Enable `LIBIGL_PREDICATES` and `LIBIGL_BUILD_TUTORIALS` in a configured libigl
build, then build `1005_SignedDistanceProfile`. It requires no viewer or OpenGL.
The normal tutorial build obtains `armadillo.obj` through `igl::tutorial_data`.

```sh
cmake --build build --target 1005_SignedDistanceProfile -j 4
build/bin/1005_SignedDistanceProfile
# Optional arguments: mesh, point count, repetitions, thread count.
build/bin/1005_SignedDistanceProfile /path/to/armadillo.obj 10000 5 4
```

The default thread count is **one**, explicitly set before any libigl call;
the optional final argument overrides it. This also overrides `IGL_NUM_THREADS`
for this executable. Run separate processes to compare thread counts because
`igl::default_num_threads` caches its first value. Eigen's own thread count is
set to one to avoid nested parallel work.

To reuse existing dependency checkouts without configuring the other tutorials,
the standalone experiment harness builds this same source:

```sh
cmake -S experiments/exact-points -B build-exact-points -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPREDICATES_SOURCE_DIR="$PREDICATES_SOURCE" \
  -DEIGEN_INCLUDE_DIR="$EIGEN_SOURCE" \
  -DCATCH2_INCLUDE_DIR="$CATCH2_SOURCE/single_include" \
  -DTUTORIAL_DATA_DIR="$TUTORIAL_DATA"
cmake --build build-exact-points --target 1005_SignedDistanceProfile -j 4
build-exact-points/1005_SignedDistanceProfile "$TUTORIAL_DATA/armadillo.obj" 10000 5 4
```

## Timing and parallelism

The output is CSV with `#` metadata/comment lines. It separates:

- **preprocess**: the sum of separately timed setup components, measured once.
  These include the AABB tree and normals/winding BVH, or the Eytzinger tree,
  array triangle adjacency and optional exposure cache.
- **query**: median/min/max over the requested repetitions, reusing all setup.
  Each of the five paths uses the **same** `parallel_queries` helper and
  `igl::parallel_for` threshold of 10,000 points, with the same thread count,
  scheduling and exception handling. It calls the scalar query kernels, which
  do not launch nested workers. The unsigned row uses the AABB query underlying
  `point_mesh_squared_distance`; the winding row retains the exact expression
  used by `signed_distance`. All normal API outputs are computed and retained.
- **total**: the actual convenience APIs for unsigned/pseudonormal/winding
  distance, and tree+adjacency(+cache) construction followed by the exact-sign
  batch API. Setup and destruction are included. All five batch paths use
  `igl::parallel_for` with the same 10,000-point threshold.

Below the threshold, or with one thread, all five query methods run serially.
The example reports whether the prepared methods actually used a thread pool
and rejects inconsistent behavior. `IGL_PARALLEL_FOR_FORCE_SERIAL`, if set at
compile time, is reported explicitly. Full-batch methods are not wrapped in an
additional parallel loop.

Every method gets an untimed warm-up. Method order rotates between repetitions.
Checksums and accuracy comparisons are evaluated outside timing. Loading,
sampling and separate mesh validation are excluded; validation intrinsic to
exposure preprocessing remains part of its setup cost. Prepared output storage
is reused, and total-call output storage is reused after the warm-up. Different
APIs return different amounts of metadata, which remains part of their costs.

“Uncached” means **no edge/vertex exposure cache**. Its AABB tree and connectivity
are still reused for the query phase. Both exact-sign paths perform one tree
traversal per query and reuse the exact selected feature for signing.

## Recorded Armadillo run

GCC 13.3, Release `-O3 -DNDEBUG`, strict FP, AMD EPYC 9B45. One thread was pinned
to CPU 24; four threads were pinned to CPUs 24–27. Five repetitions, no concurrent
compilation. Raw output is in
`experiments/exact-points/results/armadillo-{1-thread,4-threads}.csv`.

Median **query time for all 10,000 points**, with preprocessing already done:

| Method | 1 thread | 4 threads |
| --- | ---: | ---: |
| Unsigned squared distance | 87.82 ms | 22.47 ms |
| Pseudonormal signed distance | 87.58 ms | 22.28 ms |
| Fast winding signed distance | 93.80 ms | 24.04 ms |
| Exact sign, uncached | 178.89 ms | 45.76 ms |
| Exact sign, cached | 154.56 ms | 39.48 ms |

Median **full-call time including setup**:

| Method | 1 thread | 4 threads |
| --- | ---: | ---: |
| Unsigned squared distance | 139.47 ms | 77.16 ms |
| Pseudonormal signed distance | 182.17 ms | 116.63 ms |
| Fast winding signed distance | 185.90 ms | 98.08 ms |
| Exact sign, uncached | 255.84 ms | 118.54 ms |
| Exact sign, cached | 401.53 ms | 275.07 ms |

Caching reduces query time by about 14% here, but computing exposure costs
about 166 ms and is serial. It therefore increases the cost of a single batch;
its benefit depends on reuse. The observed setup/query difference amortizes
at roughly 68,000 points with one thread or 266,000 with four, for this sampling
distribution and machine. These are estimates from component timings, not
universal break-even counts. Small unsigned/pseudonormal timing differences
are within measurement variability.

All sampled sign classifications agreed. Cached and uncached exact outputs
were identical, as were each prepared path and its corresponding convenience
API. Unsigned and pseudonormal magnitudes differed from the exact-sign numerical
magnitudes by at most `1.67e-16`. The fast-winding path differed by at most
`9.04e-4`: libigl returns `sqrt(squared_distance) * (1 - 2*abs(w))`, which scales
the magnitude when the approximate winding number differs from zero or one.
The example preserves that behavior instead of replacing it with a thresholded
sign. These sampled comparisons are not a proof of accuracy for arbitrary
queries, especially near the surface.

The sample contains 600 inside and 9,400 outside points, with closest features
at 6,110 vertices, 2,982 open edges and 908 open faces. The exact-sign methods
retain their documented arithmetic range and floating-point assumptions.

Release tests passed 55,814 assertions in 21 cases, including static linkage,
serial/parallel threshold checks and worker-exception propagation. The suite
also passed with a 1 MiB stack. AddressSanitizer/UndefinedBehaviorSanitizer with
leak detection passed all three CTest targets and the full 10,000-point
Armadillo example with four threads. See
`experiments/exact-points/results/validation-armadillo-profile.json` for the
mesh revision/hash and validation metadata.

The [search-cost investigation](../../experiments/exact-points/SEARCH_PROFILE.md)
separates ordinary AABB/Eytzinger traversal, exact leaf comparisons, interval
pruning bounds and shared-feature ties using the same queries and thread policy.
