# Exact closest-point and local mesh-sign predicates

This extends `experiment/exact-rays-shewchuk` with five operations. All public
functions have `.h/.cpp` pairs under `include/igl/predicates/`; link
`igl::predicates`. The ray APIs and core ray-box implementation are unchanged.

| Function | Result |
| --- | --- |
| `exact_point_triangle_closest_feature(q,a,b,c,id)` | `ClosestPointFeature::{VERTEX,EDGE,FACE}` and local feature id |
| `exact_point_triangle_squared_distance_compare(q,a,b,c,d,e,f)` | -1 / 0 / +1 for first closer / equal / farther |
| `exact_point_mesh_closest_triangle(q,V,F,B1,B2,leaf)` | An exactly minimizing face id, choosing the smallest id on ties; -1 for empty tree |
| `exact_point_mesh_classify(q,V,F,face,TT,TTi)` | `PointMeshLocation::{INSIDE,ON_SURFACE,OUTSIDE}` (-1 / 0 / +1) |
| `exact_mesh_feature_exposure(V,F,TT,TTi,exposure)` | Query-independent edge and vertex exposure labels |

The cached signing overload is
`exact_point_mesh_classify(q,V,F,face,exposure)`. `MeshFeatureExposure` owns an
integer `edges` matrix (#F by 3) and `vertices` vector (#V). Their values are
`FeatureExposure::{FLAT,EXPOSED_TO_OUTSIDE,EXPOSED_TO_INSIDE,UNREACHABLE}`
(0 / +1 / -1 / 2). There are no retained expansions or exact constructed points.

Corner ids are 0, 1, 2. Edge **i joins i to (i+1)%3**, matching the array version
of `triangle_triangle_adjacency`, rather than the opposite-corner convention.
The open face has id 0. Triangle and tree distance queries support degenerate
triangles: among coincident representations, vertices precede open edges, and
then the lowest local feature id wins.

Signing and cache preprocessing require an embedded, closed, vertex-manifold
mesh with nondegenerate faces and consistent orientation out of the solid.
Cavity shells therefore point into their cavities. The supplied closest face
must be an actual global minimizer; any exact minimizer gives the same sign.
Preprocessing checks adjacency reciprocity/orientation, degenerate faces, and
single-cycle vertex links. It does not detect geometric self-intersections or
validate the global orientation. Unreferenced vertices are `UNREACHABLE`.
Invalidate exposure caches and bounding boxes after geometry/topology changes.

## Usage

```cpp
#include <igl/eytzinger_aabb.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>

Eigen::MatrixXd minimum(F.rows(),3), maximum(F.rows(),3), B1, B2;
for(Eigen::Index f=0; f<F.rows(); ++f)
{
  minimum.row(f)=maximum.row(f)=V.row(F(f,0));
  for(int i=1; i<3; ++i)
  {
    minimum.row(f)=minimum.row(f).cwiseMin(V.row(F(f,i)));
    maximum.row(f)=maximum.row(f).cwiseMax(V.row(F(f,i)));
  }
}
Eigen::VectorXi leaf;
igl::eytzinger_aabb(minimum,maximum,B1,B2,leaf);
Eigen::MatrixXi TT, TTi;
igl::triangle_triangle_adjacency(F,TT,TTi); // array overload
igl::predicates::MeshFeatureExposure exposure;
igl::predicates::exact_mesh_feature_exposure(V,F,TT,TTi,exposure);

const Eigen::RowVector3d q(0.1,0.2,0.3);
const Eigen::Index f=igl::predicates::exact_point_mesh_closest_triangle(q,V,F,B1,B2,leaf);
if(f>=0)
{
  const auto location=igl::predicates::exact_point_mesh_classify(q,V,F,f,exposure);
  // Alternatively: exact_point_mesh_classify(q,V,F,f,TT,TTi).
}
```

The caller builds the tree once and reuses it. Boxes must conservatively contain
their triangles; a box built from rounded approximate closest points is not a
replacement. Neither an approximate closest point nor an approximate distance
is returned by these APIs. Explicit static-library instantiations cover the
usual double row/column queries with `MatrixXd/MatrixXi/VectorXi` mesh/tree
storage, and float/double row/column leaf queries.

## Metric predicates and arithmetic contract

Let `r=q-a`, `e=b-a`, and `f=c-a`; all differences below are exact differences
of the original inputs, not previously rounded vector constructions.

- Vertex `a` is closest exactly when `r.e <= 0` and `r.f <= 0`.
- Open edge `ab` is closest when `0 < r.e < e.e` and
  `(e.e)(r.f) - (r.e)(e.f) <= 0`, after testing vertices.
- Otherwise the closest point is in the open face.

Squared distances are represented transiently as positive-denominator fractions:

| Feature | Numerator N | Denominator D |
| --- | --- | --- |
| Vertex | `r.r` | 1 |
| Edge | `(r.r)(e.e) - (r.e)^2` | `e.e` |
| Face | `(r.(e cross f))^2` | `(e cross f).(e cross f)` |

The ordering sign is `N1*D2 - N2*D1`, of degree at most ten. Proven common
vertices or common open edges give an exact tie directly, including reversed
edge orientation. Other comparisons recompute their own fractions. No exact
point, parameter, normal, or expansion persists across comparisons.

A cheap filter evaluates each expression and its absolute permanent. It accepts
only signs separated from zero by `256*DBL_EPSILON*permanent`. For a precise
roundoff count, assign weight 1 to an input difference, `max(w1,w2)+1` to an
addition/subtraction, and `w1+w2+1` to a product. This counts the rounding factors
in each fully expanded monomial. Dot products of differences have weight 5;
edge numerators 12; face numerators 17; face denominators 11; and the final
comparison at most 30. The positive permanent has the same bound. With unit
roundoff `u`, the error relative to the computed permanent is bounded by
`gamma_30/(1-gamma_30)`, where `gamma_n=n*u/(1-n*u)`. The implemented bound is
much larger, including its own rounding. Uncertain signs are retried with
outward-rounded intervals and finally Shewchuk expansions. A zero from the fast
filter is accepted only when both the value and permanent are zero.

The expansions use the existing dependency's `Two_Diff`, expansion sum, and
expansion scale, which use its `Two_Sum`/`Two_Product` primitives. Buffers are
local stack storage. A finite nonoverlapping binary64 expansion has at most
2098 nonzero components. Sum and scale scratch reserve their full input-length
bounds before checking the output length. The measured leaf calls allocate
nothing; traversal and one-ring collection use ordinary temporary vectors.

Inputs must be finite float or double, using the double predicates backend,
round-to-nearest, and no reassociation or floating-point contraction. As with
[Shewchuk's arithmetic](https://www.cs.cmu.edu/~quake/robust.html), the exactness
contract excludes intermediate overflow and loss through underflow. The degree
ten comparison can reach these limits sooner than orient3d. This is **not** an
arbitrary-exponent exact-real backend. For example, coordinates that are zero or
have absolute value between `2^-40` and `2^40` are a conservative sufficient
range for these expressions (not a necessary restriction). Tests also include
specific safe expressions outside that range. Nonfinite input coordinates throw;
not every unsupported intermediate exponent is detected dynamically.

## Conservative tree traversal

Each node gets a lower bound on squared point-box distance. Subtractions,
squares, and positive accumulation have at most five relative rounding factors.
Subtracting `16*DBL_EPSILON` times the computed sum, plus four minimum subnormal
ulps, bounds downward including rounding of the correction itself. Nonfinite
computed sums give the safe lower bound zero. Empty padding nodes do not read
uninitialized box coordinates.

The current exact minimizing triangle supplies an upper bound: the minimum of
outward-rounded squared distances to its vertices and, when its denominator
interval is positive, an outward-rounded interval quotient for its classified
feature. A node is pruned only when its lower bound is **strictly greater** than
that upper bound. Leaf candidates are compared by the metric predicate; an equal
distance still permits replacement by a smaller face id. Consequently, rounded
quotients or rounded constructed closest points never decide the winner.

## Local sign and query-free exposure

For an open face, Shewchuk `orient3d(a,b,c,q)>0` means inside for an outward
oriented triangle. A zero squared distance is checked before any exposure lookup.

For directed edge `(v,a)`, write its oriented faces `(v,a,b)` and `(v,c,a)`.
The fold determinant is `T=((a-v) cross (b-v)).(c-v)`; Shewchuk's orientation
returns the opposite sign. Thus `orient3d(v,a,b,c)>0` is convex / outside,
negative is reflex / inside, and zero uses an incident face's orientation on q.

At a closest vertex, put `d=q-v`. Every neighbor direction `e_i` satisfies
`d.e_i<=0`. Select any edge minimizing `(d.e_i)^2/(e_i.e_i)` using exact metric
comparisons. To see why ties are harmless, let u be its unit direction and
`m=d.u`, the maximum normalized dot product (all are nonpositive). For every
unit neighbor w, the supporting normal `h=d-m*u` satisfies

```
h.w <= m*(1-u.w) <= 0,       h.d = |h|^2 > 0.
```

Nondegenerate incident triangles exclude `h=0`. In an edge-normal cross-section,
both incident face rays lie away from q; q occupies the larger sector. The fold
therefore determines its side. If the edge is flat, the face plane is the
supporting plane and the query's face orientation is nonzero. The segment to
the closest vertex cannot cross the surface. The algorithm uses no unit vectors,
support-plane construction, ray, first-hit search, or symbolic perturbation.

For preprocessing, first test whether all neighbors lie in the plane of an
incident face. Such a vertex is `FLAT`. For a nonplanar ring, enumerate pairs of
neighbor directions and use orient3d to find a plane through them and v with all
neighbors in one closed half-space. All-zero pairs are collinear and skipped.
This finds a supporting hull facet whenever v is on the hull boundary, including
when v is in the **interior of a hull edge or facet**. It costs O(valence^3),
avoiding a hull dependency. No supporting plane means v is strictly inside the
3D hull and therefore `UNREACHABLE`.

A neighbor on a supporting plane supplies an incident edge and hence its fold
label. If that edge is flat, an off-plane neighbor determines the inward hull
side; the opposite side supplies the cached label. The nonzero admissible
normal cone is connected for a nonplanar ring: it is convex and contains no
opposite directions because the neighbors span 3D. Its directions cannot cross
an incident face. Therefore a single inside/outside label applies to every
noncoincident query that can have this vertex as its closest point.

## Validation and timings

The standalone harness includes the normal libigl unit test files and existing
exact-ray leaf regression tests, plus separate static-link and allocation/path
checks. Coverage includes:

- 4,000 bounded-integer cases checked against independently enumerated rational
  vertex/edge/face projections; repeated at powers-of-two scales through ±60.
- 1,000 broad-exponent and one-ulp-neighbor cases checking filter certificates
  against full expansion evaluation.
- All feature types, exact boundaries, degenerate triangles, float promotion,
  thin cancellation, and nonzero differences lost by rounded squared distances.
- 300 tree searches against exhaustive exact leaf comparisons, empty trees,
  deterministic ties, and a distant-query pruning regression.
- Convex, reflex, flat, unreachable saddle, hull-edge and hull-facet vertices;
  flat supporting edges in nonflat rings; tied closest triangles; affine
  transforms, reflections with corrected orientation, and cavity shells.
- 2,000 additional queries on 20 irregular closed stars compared with an
  independent long-double solid-angle winding calculation. This is a secondary
  numerical cross-check, not an exact oracle for points arbitrarily near a face.

See `results/` for saved test summaries and timing CSVs. Leaf allocation/path
instrumentation confirms zero allocations, zero expansion sums for an ordinary
filtered query and a shared-edge tie, and actual expansion work for a distance
comparison whose rounded distances coincide. The suite also runs with a 1 MiB
stack. The original exact-ray source and ray-box implementation are unchanged.
The final Release run passed **52,760 assertions in 16 cases**, including the
existing ray leaf tests. All three CTest targets passed in both Release and
AddressSanitizer/UndefinedBehaviorSanitizer builds, with leak detection enabled.

GCC 13.3, AMD EPYC 9B45, Release `-O3 -DNDEBUG`, strict FP, generic x86-64 code,
single thread pinned to CPU 24. The mesh benchmark uses normalized subdivided
icosahedra and 256 seeded queries uniform in `[-1.5,1.5]^3`. Each row is the
median of seven samples of at least 20 ms, with checksums and no concurrent
compilation. Tree construction and adjacency/cache preprocessing are excluded
from queries. Signing receives an already computed exact closest triangle.
The approximate baseline is `eytzinger_aabb_sdf<true>` with
`point_simplex_squared_distance<3>`; it computes an approximate distance whereas
the exact query returns a certified minimizing index, so it is a practical cost
comparison between different accuracy contracts.

| Operation | 320 faces | 5,120 faces |
| --- | ---: | ---: |
| Approximate tree query | 0.87 µs | 4.26 µs |
| Exact tree query | 2.97 µs | 8.44 µs |
| Closest-feature predicate | 31 ns | 32 ns |
| Compare nearest/adjacent triangle distances | 64 ns | 65 ns |
| Local signing, including classification | 420 ns | 494 ns |
| Cached signing, including classification | 35 ns | 36 ns |
| Exposure preprocessing, whole mesh | 0.40 ms | 6.44 ms |

The first interval-only implementation took about 1.3 µs per classification and
46/107 µs per exact tree query. `interval-only.patch` recreates that arithmetic
path and conservative box evaluation for an ablation; its saved measurements
are in `results/timings-interval-only.csv`. The cheap certified filter, direct
common-feature ties, and bounded-error box lower bound remove most of that cost.
The remaining exact tree overhead is about 3.4× / 2.0× versus approximate
traversal on these workloads. Boundary-heavy inputs needing expansions and very
high-valence preprocessing can be slower; these are not universal mesh timings.

## Reproduce

Reuse existing dependency checkouts; do not duplicate installed dependencies.
The predicates override must contain this branch's expansion API patch.

```sh
cmake -S experiments/exact-points -B build-exact-points -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPREDICATES_SOURCE_DIR="$PREDICATES_SOURCE" \
  -DEIGEN_INCLUDE_DIR="$EIGEN_SOURCE" \
  -DCATCH2_INCLUDE_DIR="$CATCH2_SOURCE/single_include"
cmake --build build-exact-points -j 4
ctest --test-dir build-exact-points --output-on-failure
(ulimit -s 1024; build-exact-points/test_exact_points)
taskset -c 24 build-exact-points/benchmark_exact_points
```

For ASan/UBSan, configure a separate build with
`-DCMAKE_CXX_FLAGS='-g -fsanitize=address,undefined -fno-omit-frame-pointer'`,
the corresponding `CMAKE_C_FLAGS`, and
`-DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'`. Build the
`test_exact_points`, `static_smoke`, and `audit_allocations` targets and run CTest.
For the interval-only ablation, apply `interval-only.patch` in a disposable
checkout, build the benchmark there, and retain strict floating-point flags.
