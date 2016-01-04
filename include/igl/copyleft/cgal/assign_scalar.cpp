// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "assign_scalar.h"

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const typename CGAL::Epeck::FT & cgal,
  CGAL::Epeck::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const typename CGAL::Epeck::FT & cgal,
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
  const double & c,
  double & d)
{
  d = c;
}
