// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/barycentric_interpolation.h>
#include <igl/readOBJ.h>
#include <igl/blue_noise.h>
#include <igl/knn.h>
#include <igl/octree.h>
#include <igl/slice.h>

//TEST_CASE("blue_noise: decimated-knight", "[igl]")
//{
//  Eigen::MatrixXd V;
//  Eigen::MatrixXi F;
//  igl::readOBJ(test_common::data_path("decimated-knight.obj"),V,F);
//  const double r = 0.01;
//  Eigen::MatrixXd B,P;
//  {
//    Eigen::VectorXi I;
//    igl::blue_noise(V,F,r,B,I,P);
//  }
//  // There should be ~4000 samples on this model
//  REQUIRE(P.rows() > 3000);
//  std::vector<std::vector<int> > point_indices;
//  Eigen::MatrixXi CH;
//  Eigen::MatrixXd CN;
//  Eigen::VectorXd W;
//  igl::octree(P,point_indices,CH,CN,W);
//  Eigen::MatrixXi I;
//  igl::knn(P,2,point_indices,CH,CN,W,I);
//  Eigen::MatrixXd P2;
//  igl::slice(P,I.col(1).eval(),1,P2);
//  Eigen::VectorXd D = (P-P2).rowwise().norm();
//  REQUIRE(D.minCoeff() > r);
//}
//
//namespace blue_noise
//{
//  template <typename DerivedA, typename DerivedB>
//  void assert_neq_different_sizes(
//    const Eigen::MatrixBase<DerivedA> & A,
//    const Eigen::MatrixBase<DerivedB> & B)
//  {
//    // test_common::assert_neq requires same sizes
//    if (A.rows() == B.rows() &&
//        A.cols() == B.cols())
//    {
//      test_common::assert_neq(A, B);
//    }
//  }
//
//  template<typename URBG>
//  void test_reproduce()
//  {
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//    igl::readOBJ(test_common::data_path("decimated-knight.obj"),V,F);
//
//    static constexpr double r = 0.1;
//    static constexpr int seed = 0;
//    Eigen::MatrixXd P1, Px1;
//    {
//      Eigen::MatrixXd B;
//      Eigen::VectorXi I;
//      URBG rng1(seed);
//      igl::blue_noise(V,F,r,B,I, P1,rng1);
//      igl::blue_noise(V,F,r,B,I,Px1,rng1);
//    }
//    Eigen::MatrixXd P2, Px2;
//    {
//      Eigen::MatrixXd B;
//      Eigen::VectorXi I;
//      URBG rng2(seed);
//      igl::blue_noise(V,F,r,B,I, P2,rng2);
//      igl::blue_noise(V,F,r,B,I,Px2,rng2);
//    }
//
//    test_common::assert_eq(P1, P2);
//    test_common::assert_eq(Px1, Px2);
//
//    assert_neq_different_sizes(P1, Px1);
//    assert_neq_different_sizes(P2, Px2);
//  }
//}
//
//
//TEST_CASE("blue_noise: minstd_rand0_reproduce", "[igl]")
//{
//  blue_noise::test_reproduce<std::minstd_rand0>();
//}
//TEST_CASE("blue_noise: minstd_rand_reproduce", "[igl]")
//{
//  blue_noise::test_reproduce<std::minstd_rand>();
//}
//TEST_CASE("blue_noise: mt19937_reproduce", "[igl]")
//{
//  blue_noise::test_reproduce<std::mt19937>();
//}
//TEST_CASE("blue_noise: mt19937_64_reproduce", "[igl]")
//{
//  blue_noise::test_reproduce<std::mt19937_64>();
//}
