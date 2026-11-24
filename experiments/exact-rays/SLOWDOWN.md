# Classification slowdown investigation

The classification regressions from the initial Shewchuk experiment are removed
on the tested GCC configurations. The mitigation preserves the public API,
thread-safe initialization, floating-point filters, arithmetic assumptions,
Shewchuk arithmetic, and local stack storage. No expansion caching was added.

## What caused the regressions

The ordinary path had a code-generation problem. Before this change, GCC emitted
the public adapter as a 192-byte frame that copied five input vectors and called
an out-of-line classifier with a 1792-byte frame. That classifier mixed cheap
sign filters with expansion scratch and external arithmetic calls. Even inputs
whose signs were immediately certified paid for this structure. The original
PR's classifier was inlined into its public template.

Separating the filtered and adaptive paths lets GCC inline the filtered work
into the public template. In the measured native build, the adapter becomes a
tail jump and the public template has a 256-byte frame. The 1792-byte adaptive
frame is entered only when a filter cannot certify a sign. The frame sizes are
evidence of changed generated code, not a claim that reserving stack space alone
accounts for the entire timing difference.

The axis-aligned coplanar path also did unnecessary arithmetic. It formed all
six determinant monomials even when each had an exactly zero factor. Skipping
those terms, and recognizing equal coordinates before calling `Two_Diff`, cuts
public C calls for the measured coplanar classification:

| Public call | Before | After |
| --- | ---: | ---: |
| `predicates_two_diff` | 6 | 2 |
| `predicates_expansion_scale` | 12 | 0 |
| `predicates_expansion_sum` | 5 | 0 |
| `orient3d` | 1 | 1 |

These are measured external calls. They do not count arithmetic inside
`orient3d`. Counts come from a separate executable using GNU ld `--wrap`, so the
instrumentation does not affect the timing executable.

The initialization guard was not a material cause. Removing it after warmup
changed random classification from 14.36 to 14.29 ns in the isolated probe.
Moving initialization into fallback branches was slower in that probe. The
production implementation retains the existing initialization behavior.

## Controlled probes and selected changes

The initial native ablations measured the following median nanoseconds per
classification. Every row is an isolated change from `4770e3f1` unless stated
otherwise; effects should not be added together.

| Probe | Random | All-hit | Axis-aligned coplanar |
| --- | ---: | ---: | ---: |
| Unmodified experiment | 14.36 | 17.31 | 163.89 |
| Remove warmed initialization guard | 14.29 | 17.36 | 163.87 |
| Bypass copies for contiguous double inputs | 13.18 | 16.37 | 163.70 |
| Outline only the determinant fallback | 13.60 | 16.78 | 168.65 |
| Skip zero monomials | 14.31 | 17.37 | 113.57 |
| Skip differences of equal coordinates | 14.33 | 17.45 | 157.17 |
| Reject behind-origin intersections before edges | 11.39 | 17.38 | 164.79 |
| Separate the entire filtered classifier (second run) | 10.42 | 11.85 | 166.21 |

The selected implementation combines:

1. A filtered classifier containing only certified nonzero signs and an
   adaptive fallback. An uncertain denominator, numerator, or edge sign retries
   through the adaptive path. That path can repeat a few cheap filters, keeping
   the common path small without complicating scratch ownership.
2. An early rejection when the exact numerator and denominator have opposite
   signs. The triangle is behind the origin, so edge tests cannot make it a hit.
   A zero numerator still reaches adaptive evaluation, preserving `t = 0` hits.
3. Exact-zero factor and equal-coordinate shortcuts in determinant expansion
   evaluation. These also help hit ordering.

Production code uses ordinary `inline`, with no compiler-specific forced-inline
or no-inline attributes. The separate probe verified that this portable version
retains the improvement. The direct-input probe was not adopted: the public
template still handles float inputs and Eigen expressions through its existing
coordinate copies. Its input contract is unchanged.

## Before and after

These results link the original PR, the frozen initial experiment, and the
mitigated implementation into one executable, with identical datasets and
compiler flags. Entries are nanoseconds per operation; lower is better.

| Native workload | Original PR | Initial experiment | Mitigated |
| --- | ---: | ---: | ---: |
| Random classification | 11.85 | 14.30 | **8.12** |
| All-hit classification | 13.60 | 17.53 | **11.75** |
| Dense all-hit classification | 14.37 | 18.80 | **12.93** |
| Axis-aligned coplanar classification | 157.46 | 166.48 | **106.72** |
| Rotated coplanar classification | 592.67 | 175.39 | **135.14** |
| Thin-triangle classification | 1219.67 | 457.97 | **223.66** |
| Ordinary hit ordering | 1465.92 | 404.43 | **261.05** |
| Dense hit ordering | 3637.45 | 2188.86 | **2185.53** |
| Thin-triangle hit ordering | 1064.39 | 384.27 | **226.13** |

The dense workload applies an invertible transform to hit triangles and the
ray, avoiding systematic zero determinant entries. Its ordering time is
effectively unchanged from the initial experiment, so the zero tests did not
introduce a measurable regression in this workload. The all-hit cases also show
that the classification improvement is not solely an effect of early rejection.

The generic build shows the same direction: random classification improves from
15.32 to 10.60 ns, all-hit classification from 18.52 to 14.84 ns, and axis-aligned
coplanar classification from 170.92 to 101.67 ns. Dense all-hit classification
improves from 19.83 to 16.08 ns. All measured classification workloads beat the
original PR in both configurations.

End-to-end mesh queries retain their gains over the original PR: about 1.2x for
first-hit random soups and 1.4–1.5x for all-hit random soups. Dense layers now run
about 2.0–2.1x faster for first hit and 4.9–5.3x faster for all hits. Tree building
and output allocation remain included. These synthetic workloads are not a
claim about every mesh or processor.

## Method and validation

Environment: AMD EPYC 9B45, GCC 13.3, Eigen 5.0.1, Linux x86-64. Diagnostic runs
were pinned to logical CPU 24. Each configuration uses `-O3 -DNDEBUG
-fno-fast-math -ffp-contract=off`; native adds `-march=native`. There was no
concurrent compilation during timing. Diagnostic measurements use nine samples
of at least 20 ms per variant, shuffling variant order each round and reporting
median/min/max. The existing mesh benchmark uses seven samples of at least
50 ms, alternating original and modified implementation order.

`baseline` means PR revision `4777048f`; `previous` means the initial experiment
at `4770e3f1`; `current` means this mitigation. The Python diagnostic script
extracts the frozen implementation from Git into disposable build directories,
gives each variant distinct symbols, and checks their results before timing.
It does not modify the production implementation while generating probes.

The expanded baseline comparison checks 125,000 classifications, 75,031 valid
hit comparisons, and 1,200 mesh queries, with zero mismatches in generic and
native builds. Native unit checks pass 21,294 assertions in 12 cases, excluding
the cube-data case. New boundary coverage exercises uncertain edge signs,
vertices, `t = 0`, behind-origin rejection, and winding changes. Existing integer
oracles, thin-triangle checks, static APIs, core-only linking, and the zero-heap
leaf audit remain in the harness. The call-count audit explicitly asserts that
the zero-product shortcut is exercised.

All five CTest checks also pass under AddressSanitizer and
UndefinedBehaviorSanitizer. The native baseline comparison passes with a 1 MiB
stack. The leaf allocation audit still observes zero allocations, versus 177
for the original PR's thin classification plus ordering pair. The known-failure
reproducers retain their previous outputs, and the core ray-box files are still
unchanged.

The original range assumptions and the separate ray-box and approximate-hit
limitations still apply; these changes do not expand the correctness contract.

## Reproduce

Build the experiment using the commands in [README.md](README.md), using
`build-slowdown` for the native configuration. Then run:

```sh
uv run --no-project --python /usr/bin/python3 experiments/exact-rays/diagnose.py \
  build-slowdown --eigen "$EIGEN_SOURCE" --predicates "$PATCHED_PREDICATES_SOURCE" \
  --variants previous
taskset -c 24 build-slowdown/diagnostics/diagnose
build-slowdown/diagnostics/predicate_calls
```

Choose an allowed CPU on a different machine. For generic diagnostics, pass
`--generic` and a build directory whose libraries were also compiled without
`-march=native`. The script targets GCC-compatible compilers and GNU ld.

Omitting `--variants` reproduces the first set of ablations. The second set used
`--variants previous force value split zero_equal_early force_zero_equal_early
value_zero_equal_early split_zero_equal_early`; the portability check used
`--variants previous split_zero_equal_early splitplain_zero_equal_early`.
`no_init` and `direct` are diagnostic-only probes with warmed initialization and
contiguous double inputs; they are not production alternatives.

Raw results are in `results/diagnostics-*.csv`,
`results/classification-mitigated-*.csv`, `results/performance-mitigated-*.csv`,
`results/accuracy-mitigated-*.csv`, and `results/predicate-calls.csv`.
