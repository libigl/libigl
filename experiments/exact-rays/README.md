# Shewchuk exact-ray experiment

The report below records the initial experiment at `4770e3f1`. A subsequent
investigation removed the observed classification regressions; see
[the slowdown investigation](SLOWDOWN.md) for the implementation changes,
controlled measurements, and latest validation results.

This branch starts at libigl PR #2549 revision
`4777048f11ac9d4d26078b32248a6f96b95460f2`. It replaces the PR's allocating
expansion implementation with Shewchuk's existing arithmetic and local stack
buffers. Expansion caching between hit comparisons is deliberately unchanged:
each comparison recomputes both triangles' numerators and denominators.

## API and implementation

The public functions now live in `igl::predicates`, with matching headers and
implementation files under `include/igl/predicates/`:

- `exact_ray_triangle_intersect`: returns `RayTriangleIntersection` (`None`,
  `Hit`, or `Coplanar`).
- `exact_ray_triangle_hit_compare`: returns the sign of `t1 - t2`, assuming both
  triangles are transverse hits.
- `exact_ray_mesh_intersect` and `exact_ray_mesh_intersect_all`: retain the face
  index and approximate `Hit` overloads from the PR.

Include the new header paths and link `igl::predicates`. The original unmerged
PR's core header paths and namespace are removed. For example:

```cmake
set(LIBIGL_PREDICATES ON CACHE BOOL "" FORCE)
add_subdirectory(libigl)
target_link_libraries(my_target PRIVATE igl::predicates)
```

```cpp
#include <igl/predicates/exact_ray_triangle_intersect.h>
#include <igl/predicates/exact_ray_triangle_hit_compare.h>

auto classification = igl::predicates::exact_ray_triangle_intersect(o,d,a,b,c);
// When both triangles have classification Hit:
int order = igl::predicates::exact_ray_triangle_hit_compare(o,d,a,b,c,e,f,g);
```

The ray-box implementation and Eytzinger traversal stay byte-for-byte unchanged
in core. They have no predicates dependency. The mesh routines still build a
tree per query and retain the original traversal and approximate-hit behavior.

The point-plane numerator sign uses the existing cheap filter, followed by
Shewchuk's `orient3d` when uncertain. Determinants containing the free direction
vector need their own high-level evaluation: forming `origin + direction` would
round the ray. Their fallbacks use the dependency's `Two_Diff`,
`scale_expansion_zeroelim`, and `fast_expansion_sum_zeroelim`. Products within
that dependency use its original `Two_Product`, `Two_Sum`, and `Fast_Two_Sum`.
There is no second implementation of those arithmetic primitives in libigl.

Hit ordering evaluates the sign of `N1*D2 - N2*D1`, after making both
denominators positive. Numerator buffers have capacity 192 and denominator
buffers 96. The accumulator uses the finite-double nonoverlap bound: at most
2098 nonzero components can occupy disjoint bit positions from -1074 through
1023. Each scaled numerator adds at most 384 components. Two buffers of
`2098 + 384` doubles therefore bound each sum without a megabyte-sized stack
frame. Overflow and underflow remain outside the arithmetic contract.

GCC's native-build stack report gives 47,696 bytes for the comparison frame and
3,072 bytes for its determinant helper (these are individual frames, not total
call-stack usage). The regression comparison also passes with a 1 MiB stack.
The allocation audit measures zero heap allocations for a thin-triangle leaf
classification plus hit comparison, versus 177 in the original PR. Mesh queries
still allocate their tree, hit containers, and other ordinary storage.

## Dependency exposure

The pinned `libigl-predicates` API did not expose expansion arithmetic. The
companion branch `experiment/expansion-arithmetic`, commit
`b434435124d29a2e87d132636fac46e1090c6c2d`, adds thin C wrappers around the
existing routines. It also fixes reads one past the input length in
`fast_expansion_sum_zeroelim`; exactly-sized caller buffers make those reads
observable to ASan. A one-element input reproducer fails against the unmodified
pinned source and passes through the patched API.

For a self-contained experiment, libigl's FetchContent recipe applies
`cmake/recipes/external/predicates-expansion.patch` to the original pinned
dependency. The patch script accepts an already-patched checkout. An explicit
`FETCHCONTENT_SOURCE_DIR_PREDICATES` override must point to a patched checkout,
because FetchContent skips patch steps for source overrides. The companion
commit is local; no published dependency revision is required to build this
branch using its default recipe.

Both the dependency and the CMake predicates target disable fast-math and
floating-point contraction. Inputs must be float or double, with the
double-precision predicates backend and round-to-nearest arithmetic.

## Accuracy and limitations

Observed on Linux x86-64 with GCC 13.3 and Eigen 5.0.1:

- All 22,259 assertions in 12 unit test cases passed, including the PR's five
  cases and existing orientation/incircle/insphere tests.
- An independent bounded-integer oracle checked 10,000 classifications and
  5,000 hit comparisons, including winding changes and scales through
  `2^-90` to `2^90`. Thin-triangle cancellation and a direction lost by rounded
  `origin + direction` have explicit tests.
- Side-by-side comparison against the original PR: 100,000 classifications,
  50,031 valid hit comparisons, and 1,200 first/all mesh queries, with zero
  mismatches. Both generic and native builds passed. This comparison uses
  random, coplanar, thin, and ordered cases with scales `2^-60` to `2^60`.
- All five CTest checks passed under AddressSanitizer and UndefinedBehaviorSanitizer:
  unit regression, baseline comparison, allocation audit, static APIs, and
  core-only ray-box traversal. Static checks exercise row/column double inputs,
  all mesh overloads, and float leaf APIs.

The earlier review findings are preserved in `results/known-failures.csv`:

- Distant-origin and tiny-direction rays can still be culled incorrectly by
  the unchanged core ray-box code.
- The unchanged approximate `Hit` overload can still return `t = 0` for an
  exact hit at `t = 1` when its rounded determinant vanishes.
- Extreme exponents can overflow or underflow expansion arithmetic. Both
  versions misorder the explicit `2^220` and `2^-220` examples. At scale
  `2^400`, a behind-origin triangle is rejected by the baseline but classified
  as a hit by this branch. This is outside the supported arithmetic range and
  demonstrates that this experiment does **not** preserve every result for
  arbitrary finite double inputs.

These checks establish parity on the tested domain, not a proof that the mesh
query is exact for every input. The unrelated traversal and approximate-output
bugs are not fixed here.

## Performance

Measurements are single-threaded on an AMD EPYC 9B45 in a container with a
four-core quota. Both implementations are linked into the same executable,
compiled at `-O3 -DNDEBUG -fno-fast-math -ffp-contract=off`. Native builds add
`-march=native`, allowing the baseline's explicit `std::fma` to use hardware
FMA. Each row is the median of seven samples, each calibrated to at least
50 ms, alternating which implementation runs first. Checksums prevent dead-code
elimination. No build was running during the reported timing runs.

Mesh workloads are synthetic random triangle soups and dense layers, at 128
and 1024 faces, with eight rays per batch. Each timed mesh call includes tree
construction and output allocation. These results do not isolate BVH traversal
or predict performance on every production mesh.

| Workload | Generic speedup | Native speedup |
| --- | ---: | ---: |
| Random classification | 0.98x | 0.83x |
| Coplanar classification | 1.01x | 0.92x |
| Thin-triangle classification | 2.67x | 2.72x |
| Hit ordering, ordinary hits | 3.39x | 3.65x |
| Hit ordering, thin triangles | 2.68x | 2.80x |
| First hit, random meshes | 1.16–1.19x | 1.19–1.23x |
| All hits, random meshes | 1.43–1.44x | 1.50–1.52x |
| First hit, layered meshes | 1.53–1.56x | 1.52–1.56x |
| All hits, layered meshes | 2.22–2.28x | 2.15–2.19x |

Speedup is baseline time divided by experiment time. Native ordinary
classification regresses from 11.85 ns to 14.19 ns (about 20%); coplanar
classification regresses from 151.71 ns to 164.62 ns (about 9%). The initial
version called `orient3d` unconditionally and took 16.53 ns on the random
classification workload. Restoring the cheap filter reduced that overhead,
but did not eliminate the regression. Exact initialization checks, function
boundaries, and code generation remain differences from the fully inline
baseline; their individual costs have not been isolated. The zero-allocation
exact paths and mesh workloads show substantial gains. This is a performance
tradeoff, not a universal speedup.

Raw per-operation nanoseconds are in `results/performance-generic.csv` and
`results/performance-native.csv`.

## Reproduce

Use a baseline checkout at the revision above, the patched predicates checkout,
Eigen 5.0.1, and Catch2 v2's single-include directory. Reuse existing checkouts.
From this branch's root:

```sh
cmake -S experiments/exact-rays -B build-experiment -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBASELINE_SOURCE_DIR="$BASELINE_SOURCE" \
  -DPREDICATES_SOURCE_DIR="$PATCHED_PREDICATES_SOURCE" \
  -DEIGEN_INCLUDE_DIR="$EIGEN_SOURCE" \
  -DCATCH2_INCLUDE_DIR="$CATCH2_SOURCE/single_include"
cmake --build build-experiment -j 4
ctest --test-dir build-experiment --output-on-failure
build-experiment/compare_exact_rays --accuracy
build-experiment/compare_exact_rays --performance
build-experiment/compare_exact_rays --known-failures
build-experiment/audit_allocations
```

CTest omits the original cube-mesh test so it needs no test-data download. To
include that case, copy `cube.obj` from `libigl-tests-data` to
`build-experiment/data/cube.obj`, then run `build-experiment/test_exact_rays`.

For native timings, configure a separate `build-native` with
`-DCMAKE_C_FLAGS_RELEASE="-O3 -DNDEBUG -march=native"` and the same value for
`CMAKE_CXX_FLAGS_RELEASE`. Build target `compare_exact_rays`, then run its
accuracy and performance modes. For sanitizer checks, configure a separate
`build-sanitize` with both release flag variables set to
`-O1 -g -DNDEBUG -fsanitize=address,undefined -fno-omit-frame-pointer`, build all
targets, and run CTest. `-DNDEBUG` avoids the test runner's debug-mode floating
exception traps, which are incompatible with the intentional range reproducers.
