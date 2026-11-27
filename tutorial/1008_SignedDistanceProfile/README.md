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
build, then build `1008_SignedDistanceProfile`. It requires no viewer or OpenGL.
The normal tutorial build obtains `armadillo.obj` through `igl::tutorial_data`.

```sh
cmake --build build --target 1008_SignedDistanceProfile -j 4
build/bin/1008_SignedDistanceProfile
# Optional arguments: mesh, point count, repetitions, thread count.
build/bin/1008_SignedDistanceProfile /path/to/armadillo.obj 10000 5 4
```

The default thread count is **one**, explicitly set before any libigl call;
the optional final argument overrides it. This also overrides `IGL_NUM_THREADS`
for this executable. Run separate processes to compare thread counts because
`igl::default_num_threads` caches its first value. Eigen's own thread count is
set to one to avoid nested parallel work.

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

## Interpreting comparisons

Use query rows to compare methods when preprocessing can be reused. Use total
rows for a single batch on a new mesh. Exposure caching saves the local signing
step; it does not remove the closest-point search, and its setup cost can exceed
that saving for one batch. Measure both rows for the expected workload.

The example reports sign agreement, differences in numerical distance and
closest-feature counts after timing. The exact-sign API returns Eigen arrays
`S, I, bary` for signed distances, triangle indices and barycentrics. The
outputs use `Eigen::PlainObjectBase`; closest-point coordinates can be
reconstructed from the triangle vertices and barycentrics when needed. The
sign/zeros of `S` retain exact inside/on-surface/outside
classification; the zero pattern of each barycentric row retains the exact
closest feature (one positive weight for a vertex, two for an open edge, three
for an open face). Nonzero distances or active weights that round to zero are
replaced by the smallest positive value of the output scalar type to preserve
this information. Cached and uncached exact locations must agree. Fast winding
returns `sqrt(squared_distance) * (1 - 2*abs(w))`, so its
magnitude changes when the approximate winding number differs from zero or one.
These sampled comparisons do not prove accuracy for arbitrary queries.

For reproducible timings, use a Release build, keep thread counts equal, and
run several separate processes without concurrent builds or other CPU work.
The CSV includes the seed, mesh size, thread count and repetition count.
