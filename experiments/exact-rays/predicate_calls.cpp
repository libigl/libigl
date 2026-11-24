#include "workloads.h"
#include <igl/predicates/exactinit.h>
#include <predicates.h>
#include <iostream>
#include <stdexcept>

int probe_previous(const Ray &,const Triangle &);
namespace
{
struct Counts { int differences=0, scales=0, sums=0, orientations=0; } counts;
}
extern "C"
{
void __real_predicates_two_diff(double,double,double *,double *);
int __real_predicates_expansion_scale(int,const double *,double,double *);
int __real_predicates_expansion_sum(int,const double *,int,const double *,double *);
double __real_orient3d(double *,double *,double *,double *);
void __wrap_predicates_two_diff(double a,double b,double *value,double *tail)
{
  ++counts.differences;
  __real_predicates_two_diff(a,b,value,tail);
}
int __wrap_predicates_expansion_scale(int length,const double *input,double factor,double *output)
{
  ++counts.scales;
  return __real_predicates_expansion_scale(length,input,factor,output);
}
int __wrap_predicates_expansion_sum(int left_length,const double *left,int right_length,
  const double *right,double *output)
{
  ++counts.sums;
  return __real_predicates_expansion_sum(left_length,left,right_length,right,output);
}
double __wrap_orient3d(double *a,double *b,double *c,double *origin)
{
  ++counts.orientations;
  return __real_orient3d(a,b,c,origin);
}
}
int main()
{
  igl::predicates::exactinit();
  std::cout<<"workload,implementation,two_diff,expansion_scale,expansion_sum,orient3d\n";
  for(const std::string kind:{"random","ordered","coplanar","thin"})
  {
    const auto sample=exact_ray_experiment::cases(kind,1).front();
    counts={};
    const int before=probe_previous(sample.ray,sample.first);
    const Counts previous=counts;
    counts={};
    const int after=shewchuk::classify(sample.ray,sample.first);
    if(before!=after) { throw std::runtime_error("Classification changed"); }
    if(kind=="coplanar" && (previous.scales==0 || counts.scales!=0 || counts.sums!=0))
      throw std::runtime_error("Zero-product shortcut was not exercised");
    std::cout<<kind<<",previous,"<<previous.differences<<","<<previous.scales<<","
      <<previous.sums<<","<<previous.orientations<<"\n";
    std::cout<<kind<<",current,"<<counts.differences<<","<<counts.scales<<","
      <<counts.sums<<","<<counts.orientations<<"\n";
  }
}
