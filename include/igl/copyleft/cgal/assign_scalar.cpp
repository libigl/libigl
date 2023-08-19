// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "assign_scalar.h"

template <>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & rhs,
  const bool & slow_and_more_precise,
  double & lhs)
{
  if(slow_and_more_precise)
  {
    return assign_scalar(rhs,lhs);
  }else
  {
    // While this is significantly faster (100x), this does not guarantee that
    // two equivalent rationals produce the same double (e.g.,
    // CGAL::to_double(⅓) ≠ CGAL::to_double(1-⅔))
    // https://github.com/CGAL/cgal/discussions/6000 For remesh_intersections,
    // `unique` is called _after_ rounding to floats, avoiding more expensive
    // rational equality tests. To operate correctly, we need that a=b ⇒
    // double(a)=double(b). Alternatively, we could require that
    // remesh_intersections conduct its `unique` operation on rationals. This is
    // even more expensive (4x) and most of the time probably overkill, though
    // it is argueably more correct in terms of producing the correct topology
    // (despite degeneracies that appear during rounding). On the other hand,
    // this rounding-before-unique only occurs if the requested output is float,
    // so its a question of combinatorial vs geometric degeneracies. Argueably,
    // combinatorial degeneracies are more reliably detected.
    //
    //lhs = CGAL::to_double(rhs);
    lhs = CGAL::to_double(rhs.exact());
  }
}

template <>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & rhs,
  const bool & slow_and_more_precise,
  float & lhs)
{
  double d;
  igl::copyleft::cgal::assign_scalar(rhs,slow_and_more_precise,d);
  lhs = float(d);
}

// If we haven't specialized the types then `slow_and_more_precise` doesn't make sense.
template <typename RHS, typename LHS>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const RHS & rhs,
  const bool & slow_and_more_precise,
  LHS & lhs)
{
  return assign_scalar(rhs,lhs);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & cgal,
  CGAL::Epeck::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & _cgal,
  double & d)
{
  // FORCE evaluation of the exact type otherwise interval might be huge.
  const CGAL::Epeck::FT cgal = _cgal.exact();
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & _cgal,
  float& d)
{
  // FORCE evaluation of the exact type otherwise interval might be huge.
  const CGAL::Epeck::FT cgal = _cgal.exact();
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const double & c,
  double & d)
{
  d = c;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const float& c,
  float& d)
{
  d = c;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const float& c,
  double& d)
{
  d = c;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  double & d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  float& d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

#ifndef WIN32

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  CGAL::Simple_cartesian<mpq_class>::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  double & d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  float& d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

#endif // WIN32

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::copyleft::cgal::assign_scalar<float, double>(float const&, bool const&, double&);
// generated by autoexplicit.sh
template void igl::copyleft::cgal::assign_scalar<float, CGAL::Epeck::FT >(float const&, bool const&, CGAL::Epeck::FT&);
template void igl::copyleft::cgal::assign_scalar<double, double>(double const&, bool const&, double&);
template void igl::copyleft::cgal::assign_scalar<double, CGAL::Epeck::FT >(double const&, bool const&, CGAL::Epeck::FT&);
template void igl::copyleft::cgal::assign_scalar<CGAL::Epeck::FT, CGAL::Epeck::FT >(CGAL::Epeck::FT const&, bool const&, CGAL::Epeck::FT&);
#endif
