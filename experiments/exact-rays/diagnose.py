"""Build isolated classification ablations from the original experiment commit."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("build", type=Path)
parser.add_argument("--eigen", required=True)
parser.add_argument("--predicates", required=True)
parser.add_argument("--generic", action="store_true")
parser.add_argument("--variants", nargs="+")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build.resolve()
scratch = build / "diagnostics"
scratch.mkdir(exist_ok=True)


def original(name):
    return subprocess.check_output([
        "git", "-C", str(root), "show",
        "4770e3f1:include/igl/predicates/" + name], text=True, timeout=30)


header = original("exact_ray_triangle_intersect.h").replace(
    '"../igl_inline.h"', '<igl/igl_inline.h>')
implementation = original("exact_ray_triangle_intersect.cpp").replace(
    '"exactinit.h"', '<igl/predicates/exactinit.h>')
kernel = original("exact_ray_triangle_kernel.h")


def transform(name):
    code, api = kernel, implementation
    if "no_init" in name or "lazy" in name:
        api = api.replace("  igl::predicates::exactinit();\n", "")
    if "direct" in name:
        start = api.index("  const double O[3]")
        end = api.index("  return static_cast", start)
        api = api[:start] + """  const double *O = source.derived().data();
  const double *dd = dir.derived().data();
  const double *a = A.derived().data();
  const double *b = B.derived().data();
  const double *c = C.derived().data();
""" + api[end:]
    if "outline" in name:
        old = """    double expansion[192];
    const int length = det3_expansion(c0, c1, c2, expansion);
    return sign(expansion[length - 1]);"""
        assert old in code
        code = code.replace(old, "    return det3_exact(c0, c1, c2);")
        helper = """  __attribute__((noinline)) inline int det3_exact(
    const Column &c0, const Column &c1, const Column &c2)
  {
""" + old + "\n  }\n\n"
        code = code.replace("  inline int det3_sign(", helper + "  inline int det3_sign(")
    if "lazy" in name:
        code = '#include <igl/predicates/exactinit.h>\n' + code
        code = code.replace("    if(filtered != 0) { return filtered; }",
            "    if(filtered != 0) { return filtered; }\n    igl::predicates::exactinit();")
        code = code.replace("      double pa[3], pb[3], pc[3], po[3];",
            "      igl::predicates::exactinit();\n      double pa[3], pb[3], pc[3], po[3];")
    if "zero" in name:
        code = code.replace("      const int *axes = permutations[permutation];", """      const int *axes = permutations[permutation];
      if((lengths[axes[0]][0] == 1 && entries[axes[0]][0][0] == 0) ||
         (lengths[axes[1]][1] == 1 && entries[axes[1]][1][0] == 0) ||
         (lengths[axes[2]][2] == 1 && entries[axes[2]][2][0] == 0))
      {
        continue;
      }""")
    if "equal" in name:
        code = code.replace("    double value, tail;", """    if(column.p[axis] == column.q[axis])
    {
      result[0] = 0;
      return 1;
    }
    double value, tail;""")
    if "early" in name:
        code = code.replace("    if(denominator == 0) { return numerator == 0 ? 2 : 0; }",
            "    if(denominator == 0) { return numerator == 0 ? 2 : 0; }\n"
            "    if(numerator == -denominator) { return 0; }")
    if "force" in name:
        code = code.replace("  inline int classify(",
            "  __attribute__((always_inline)) inline int classify(")
    if "value" in name:
        code = code.replace("""    const double *p;
    const double *q;""", """    double p[3], q[3];
    bool diff;
    Column(const double *first, const double *second): q{0,0,0}, diff(second != nullptr)
    {
      std::copy_n(first,3,p);
      if(second) { std::copy_n(second,3,q); }
    }""")
        code = code.replace("if(!column.q)", "if(!column.diff)")
        code = code.replace("col[c]->q ?", "col[c]->diff ?")
    if "split" in name:
        code = code.replace("  inline int classify(",
            "  inline int classify_fallback(" if "plain" in name else
            "  __attribute__((noinline)) inline int classify_fallback(")
        fast = """  inline int classify(const double *origin, const double *direction,
    const double *a, const double *b, const double *c)
  {
    const Column ray{direction,nullptr}, ab{b,a}, ac{c,a};
    const int denominator = det3_filter(ray,ab,ac);
    if(denominator != 0)
    {
      const int numerator = det3_filter(ab,ac,{a,origin});
      if(numerator != 0)
      {
        EARLY_REJECTION
        const int first = det3_filter(ray,{a,origin},{b,origin});
        if(first == -denominator) { return 0; }
        if(first != 0)
        {
          const int second = det3_filter(ray,{b,origin},{c,origin});
          if(second == -denominator) { return 0; }
          if(second != 0)
          {
            const int third = det3_filter(ray,{c,origin},{a,origin});
            if(third == -denominator) { return 0; }
            if(third != 0) { return numerator == denominator ? 1 : 0; }
          }
        }
      }
    }
    return classify_fallback(origin,direction,a,b,c);
  }

"""
        fast = fast.replace("EARLY_REJECTION", "if(numerator == -denominator) { return 0; }" if "early" in name else "")
        code = code.replace("  inline void hit_fraction(", fast + "  inline void hit_fraction(")
    return code, api


names = ["previous", "no_init", "direct", "outline", "lazy", "zero", "equal", "early",
         "outline_lazy", "outline_lazy_zero_equal", "outline_lazy_zero_equal_early"]
if args.variants:
    names = args.variants
if "previous" not in names:
    names.insert(0, "previous")
flags = ["-std=c++17", "-O3", "-DNDEBUG", "-fno-fast-math", "-ffp-contract=off", "-fstack-usage"]
if not args.generic:
    flags += ["-march=native"]
includes = ["-I" + str(root / "include"), "-I" + str(root / "experiments/exact-rays"),
            "-I" + args.eigen, "-I" + args.predicates]


def compile_variant(name):
    path = scratch / name
    path.mkdir(exist_ok=True)
    code, api = transform(name)
    (path / "exact_ray_triangle_kernel.h").write_text(code)
    (path / "exact_ray_triangle_intersect.h").write_text(header)
    (path / "exact_ray_triangle_intersect.cpp").write_text(api)
    compare_source = ""
    if name == "previous":
        for suffix in ("h", "cpp"):
            comparison = original("exact_ray_triangle_hit_compare." + suffix)
            comparison = comparison.replace('"../igl_inline.h"', '<igl/igl_inline.h>')
            comparison = comparison.replace('"exactinit.h"', '<igl/predicates/exactinit.h>')
            (path / ("exact_ray_triangle_hit_compare." + suffix)).write_text(comparison)
        compare_source = '''
#define exact_ray_triangle_hit_compare compare_previous
#include "exact_ray_triangle_hit_compare.h"
int probe_compare_previous(const Ray &ray,const Triangle &first,const Triangle &second)
{
  return igl::predicates::compare_previous(ray.origin,ray.direction,
    first.a,first.b,first.c,second.a,second.b,second.c);
}
'''
    source = path / "adapter.cpp"
    source.write_text(f'''#include "api.h"
#define exact_ray_triangle_kernel kernel_{name}
#define exact_ray_triangle_intersect classify_{name}
#include "exact_ray_triangle_intersect.h"
int probe_{name}(const Ray &ray, const Triangle &triangle)
{{
  return int(igl::predicates::classify_{name}(ray.origin,ray.direction,
    triangle.a,triangle.b,triangle.c));
}}
''' + compare_source)
    obj = path / "adapter.o"
    subprocess.run(["c++", *flags, *includes, "-c", str(source), "-o", str(obj)],
                   check=True, timeout=120)
    return str(obj)


with ThreadPoolExecutor(max_workers=4) as workers:
    objects = list(workers.map(compile_variant, names))

declarations = "\n".join(f"int probe_{name}(const Ray &, const Triangle &);" for name in names)
entries = ",\n".join('{"' + name + '",+[](const Ray &r,const Triangle &t,const Triangle &){return probe_'
                      + name + '(r,t);}}' for name in names)
driver = scratch / "driver.cpp"
driver.write_text('''#include "workloads.h"
#include <igl/predicates/exactinit.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
''' + declarations + '''
int probe_compare_previous(const Ray &,const Triangle &,const Triangle &);
using Operation = int (*)(const Ray &,const Triangle &,const Triangle &);
struct Variant { const char *name; Operation operation; int repetitions=1; std::vector<double> times; };
volatile int checksum;
double duration(Operation operation,const std::vector<exact_ray_experiment::Case> &samples,int repetitions)
{
  int sum=0;
  const auto start=std::chrono::steady_clock::now();
  for(int i=0;i<repetitions;++i)
    for(const auto &sample:samples) sum+=operation(sample.ray,sample.first,sample.second);
  checksum=sum;
  return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
int main()
{
  // no_init is an attribution probe only; all timings start with initialized C state.
  igl::predicates::exactinit();
  std::cout<<"workload,variant,median_ns,min_ns,max_ns\\n"<<std::setprecision(10);
  for(const std::string kind:{"random","ordered","dense_ordered","coplanar","rotated_coplanar","thin",
                             "compare_ordered","compare_dense_ordered","compare_thin"})
  {
    const bool comparison=kind.rfind("compare_",0)==0;
    auto samples=exact_ray_experiment::cases(kind=="rotated_coplanar"?"coplanar":
      (comparison?kind.substr(8):kind),1024);
    if(kind=="rotated_coplanar")
      for(auto &sample:samples)
      {
        sample.ray.origin.z()=sample.ray.origin.x()+sample.ray.origin.y();
        sample.ray.direction.z()=sample.ray.direction.x()+sample.ray.direction.y();
        for(Point *point:{&sample.first.a,&sample.first.b,&sample.first.c}) point->z()=point->x()+point->y();
      }
    std::vector<Variant> variants={{"baseline",+[](const Ray &r,const Triangle &t,const Triangle &){return baseline::classify(r,t);}},
      {"current",+[](const Ray &r,const Triangle &t,const Triangle &){return shewchuk::classify(r,t);}},
''' + entries + '''};
    if(comparison) variants={{"baseline",baseline::compare},{"current",shewchuk::compare},
      {"previous",probe_compare_previous}};
    for(auto &variant:variants)
    {
      for(const auto &sample:samples)
        if(variant.operation(sample.ray,sample.first,sample.second)!=(comparison?
            baseline::compare(sample.ray,sample.first,sample.second):baseline::classify(sample.ray,sample.first)))
          throw std::runtime_error(std::string("Classification mismatch: ")+variant.name);
      while(duration(variant.operation,samples,variant.repetitions)<0.02) variant.repetitions*=2;
    }
    std::vector<int> order(variants.size()); std::iota(order.begin(),order.end(),0);
    std::mt19937 generator(719);
    for(int repetition=0;repetition<9;++repetition)
    {
      std::shuffle(order.begin(),order.end(),generator);
      for(int index:order)
      {
        auto &variant=variants[index];
        variant.times.push_back(duration(variant.operation,samples,variant.repetitions)*1e9
          /(variant.repetitions*samples.size()));
      }
    }
    for(auto &variant:variants)
    {
      std::sort(variant.times.begin(),variant.times.end());
      std::cout<<kind<<","<<variant.name<<","<<variant.times[4]<<","<<variant.times.front()
        <<","<<variant.times.back()<<"\\n"<<std::flush;
    }
  }
}
''')
subprocess.run(["c++", *flags, *includes, str(driver), *objects,
                str(build / "libbaseline_adapter.a"), str(build / "libshewchuk_adapter.a"),
                str(build / "lib/libpredicates.a"), "-o", str(scratch / "diagnose")],
               check=True, timeout=120)
print(scratch / "diagnose")
if "previous" in names:
    wrapped = ["predicates_two_diff", "predicates_expansion_scale", "predicates_expansion_sum", "orient3d"]
    subprocess.run(["c++", *flags, *includes,
                    str(root / "experiments/exact-rays/predicate_calls.cpp"),
                    str(scratch / "previous/adapter.o"), str(build / "libshewchuk_adapter.a"),
                    str(build / "lib/libpredicates.a"), *["-Wl,--wrap=" + name for name in wrapped],
                    "-o", str(scratch / "predicate_calls")], check=True, timeout=120)
