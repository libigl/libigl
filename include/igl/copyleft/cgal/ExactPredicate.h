// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Qingan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_COPYLEFT_CGAL_EXACT_PREDICATE
#define IGL_COPYLEFT_CGAL_EXACT_PREDICATE

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      template<typename Scalar>
      class ExactPredicate {
        public:
          typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
          typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
          typedef typename std::conditional<std::is_same<Scalar, Epeck::FT>::value,
                  Epeck, Epick>::type Kernel;

          static short orientation(const Scalar pa[2], const Scalar pb[2], const Scalar pc[2])
          {
            switch(CGAL::orientation(
                  typename Kernel::Point_2(pa[0], pa[1]),
                  typename Kernel::Point_2(pb[0], pb[1]),
                  typename Kernel::Point_2(pc[0], pc[1]))) {
              case CGAL::LEFT_TURN:
                return 1;
              case CGAL::RIGHT_TURN:
                return -1;
              case CGAL::COLLINEAR:
                return 0;
              default:
                throw "Invalid orientation";
            }
          }

          static short orientation(const Scalar pa[3], const Scalar pb[3], const Scalar pc[3], const Scalar pd[3])
          {
            switch(CGAL::orientation(
                  typename Kernel::Point_3(pa[0], pa[1], pa[2]),
                  typename Kernel::Point_3(pb[0], pb[1], pb[2]),
                  typename Kernel::Point_3(pc[0], pc[1], pc[2]),
                  typename Kernel::Point_3(pd[0], pd[1], pd[2]))) {
              case CGAL::POSITIVE:
                return 1;
              case CGAL::NEGATIVE:
                return -1;
              case CGAL::COPLANAR:
                return 0;
              default:
                throw "Invalid orientation";
            }
          }

          static short incircle(const Scalar pa[2], const Scalar pb[2], const Scalar pc[2], const Scalar pd[2])
          {
            switch(CGAL::side_of_oriented_circle(
                  typename Kernel::Point_2(pa[0], pa[1]),
                  typename Kernel::Point_2(pb[0], pb[1]),
                  typename Kernel::Point_2(pc[0], pc[1]),
                  typename Kernel::Point_2(pd[0], pd[1]))) {
              case CGAL::ON_POSITIVE_SIDE:
                return 1;
              case CGAL::ON_NEGATIVE_SIDE:
                return -1;
              case CGAL::ON_ORIENTED_BOUNDARY:
                return 0;
              default:
                throw "Invalid incircle result";
            }
          }

          static short insphere(const Scalar pa[3], const Scalar pb[3], const Scalar pc[3], const Scalar pd[3],
              const Scalar pe[3])
          {
            switch(CGAL::side_of_oriented_sphere(
                  typename Kernel::Point_3(pa[0], pa[1], pa[2]),
                  typename Kernel::Point_3(pb[0], pb[1], pb[2]),
                  typename Kernel::Point_3(pc[0], pc[1], pc[2]),
                  typename Kernel::Point_3(pd[0], pd[1], pd[2]),
                  typename Kernel::Point_3(pe[0], pe[1], pe[2]))) {
              case CGAL::ON_POSITIVE_SIDE:
                return 1;
              case CGAL::ON_NEGATIVE_SIDE:
                return -1;
              case CGAL::ON_ORIENTED_BOUNDARY:
                return 0;
              default:
                throw "Invalid incircle result";
            }
          }
      };
    }
  }
}

#endif
