// Benchmark-only full distance enclosures; production needs only upper bounds.
#pragma once
#include <igl/predicates/exact_point_triangle_kernel.h>
namespace filter_profile
{
namespace k=igl::predicates::exact_point_triangle_kernel;
using Fraction=k::Fraction<k::Filter>;
inline bool upper_bound(const Fraction &fraction,double &upper)
{ return k::distance_upper_bound(fraction,upper); }
inline int feature_distance(const k::Point &q,const std::array<k::Point,3> &p,
  int &id,Fraction &distance)
{ return k::filtered_feature_distance(q,p,id,distance); }
inline bool distance_bounds(const Fraction &fraction,k::Interval &bounds)
{
  double numerator_error,denominator_error;
  if(!k::filter_error(fraction.numerator,numerator_error) || !k::filter_error(fraction.denominator,denominator_error))
    return false;
  const double denominator_lo=k::down(fraction.denominator.value-denominator_error);
  const double denominator_hi=k::up(fraction.denominator.value+denominator_error);
  if(denominator_lo<=0) return false;
  const double numerator_lo=std::max(0.0,k::down(fraction.numerator.value-numerator_error));
  const double numerator_hi=k::up(fraction.numerator.value+numerator_error);
  bounds={std::max(0.0,k::down(numerator_lo/denominator_hi)),k::up(numerator_hi/denominator_lo)};
  return numerator_hi>=0 && std::isfinite(bounds.hi);
}

}
