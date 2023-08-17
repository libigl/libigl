// This file is part of libigl, a simple c++ geometry processing library.
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/readOBJ.h>
#include <igl/random_points_on_mesh.h>

TEST_CASE("random_points_on_mesh: decimated-knight", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(test_common::data_path("decimated-knight.obj"),V,F);
  static constexpr int n = 1000;
  Eigen::MatrixXd P;
  {
    Eigen::MatrixXd B;
    Eigen::VectorXi I;
    igl::random_points_on_mesh(n,V,F,B,I,P);
  }
  REQUIRE(P.rows() == n);
}

namespace random_points_on_mesh
{
  template<typename URBG>
  void test_reproduce()
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(test_common::data_path("decimated-knight.obj"),V,F);

    static constexpr int n = 100;
    static constexpr int seed = 0;
    Eigen::MatrixXd P1, Px1;
    {
      Eigen::MatrixXd B;
      Eigen::VectorXi I;
      URBG rng1(seed);
      igl::random_points_on_mesh(n,V,F,B,I,P1 ,rng1);
      igl::random_points_on_mesh(n,V,F,B,I,Px1,rng1);
    }
    Eigen::MatrixXd P2, Px2;
    {
      Eigen::MatrixXd B;
      Eigen::VectorXi I;
      URBG rng2(seed);
      igl::random_points_on_mesh(n,V,F,B,I,P2 ,rng2);
      igl::random_points_on_mesh(n,V,F,B,I,Px2,rng2);
    }

    test_common::assert_eq(P1, P2);
    test_common::assert_eq(Px1, Px2);

    test_common::assert_neq(P1, Px1);
    test_common::assert_neq(P2, Px2);
  }
}


TEST_CASE("random_points_on_mesh: minstd_rand0_reproduce", "[igl]")
{
  random_points_on_mesh::test_reproduce<std::minstd_rand0>();
}
TEST_CASE("random_points_on_mesh: minstd_rand_reproduce", "[igl]")
{
  random_points_on_mesh::test_reproduce<std::minstd_rand>();
}
TEST_CASE("random_points_on_mesh: mt19937_reproduce", "[igl]")
{
  random_points_on_mesh::test_reproduce<std::mt19937>();
}
TEST_CASE("random_points_on_mesh: mt19937_64_reproduce", "[igl]")
{
  random_points_on_mesh::test_reproduce<std::mt19937_64>();
}
