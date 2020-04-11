// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <test_common.h>
#include <igl/read_triangle_mesh.h>
#include <igl/winding_number.h>
#include <igl/fast_winding_number.h>
#include <igl/octree.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>

TEST_CASE("fast_winding_number: one_point_cloud", "[igl]")
{
  Eigen::MatrixXd P(1,3);
  P<<0.1,0.2,0.3;
  // Unit normal
  Eigen::MatrixXd N(1,3);
  N<<4,5,6;
  N /= N.row(0).norm();
  Eigen::VectorXd A(1,1);
  A(0) = 0.7;

  std::vector<std::vector<int > > O_PI;
  Eigen::MatrixXi O_CH;
  Eigen::MatrixXd O_CN;
  Eigen::VectorXd O_W;
  igl::octree(P,O_PI,O_CH,O_CN,O_W);

  Eigen::MatrixXd O_CM;
  Eigen::VectorXd O_R;
  Eigen::MatrixXd O_EC;
  igl::fast_winding_number(P,N,A,O_PI,O_CH,2,O_CM,O_R,O_EC);
  Eigen::MatrixXd Q(4,3);
  Q<<
     0, 0, 0,
     4,-3, 2,
    -1, 1, 3,
     0, 5,-2;
  Eigen::VectorXd WiP;
  igl::fast_winding_number(P,N,A,O_PI,O_CH,O_CM,O_R,O_EC,Q,2,WiP);
  Eigen::VectorXd WiP_cached(4);
  WiP_cached<<
    0.38779369004261133,
    -0.00041235296362485,
    -0.00362978253577090,
    -0.00041235296362485;
  test_common::assert_near(WiP,WiP_cached,1e-15);
}

TEST_CASE("fast_winding_number: meshes", "[igl]" "[slow]")
{
  const auto test_case = [](const std::string &param)
  {
    INFO(param);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path(param),V,F);
    // vertex centroid will be our query
    Eigen::MatrixXd Q = V.array().colwise().mean();

    Eigen::VectorXd Wexact(1,1);
    Wexact(0,0) = igl::winding_number(V,F,Eigen::RowVector3d(Q));

    // SOUP
    {
      INFO("soup");
      igl::FastWindingNumberBVH fwn_bvh;
      igl::fast_winding_number(V,F,2,fwn_bvh);
      Eigen::VectorXd Wfwn_soup;
      igl::fast_winding_number(fwn_bvh,2,Q,Wfwn_soup);
      test_common::assert_near(Wfwn_soup,Wexact,1e-2);
    }

    // CLOUD
    // triangle barycenters, normals and areas will be our point cloud
    {
      INFO("cloud");
      Eigen::MatrixXd BC,N;
      Eigen::VectorXd A;
      igl::barycenter(V,F,BC);
      igl::per_face_normals(V,F,N);
      igl::doublearea(V,F,A);
      A *= 0.5;
      Eigen::VectorXd Wfwn_cloud;
      std::vector<std::vector<int > > O_PI;
      Eigen::MatrixXi O_CH;
      Eigen::MatrixXd O_CN;
      Eigen::VectorXd O_W;
      igl::octree(BC,O_PI,O_CH,O_CN,O_W);
      Eigen::MatrixXd O_CM;
      Eigen::VectorXd O_R;
      Eigen::MatrixXd O_EC;
      igl::fast_winding_number(BC,N,A,O_PI,O_CH,2,O_CM,O_R,O_EC);
      igl::fast_winding_number(BC,N,A,O_PI,O_CH,O_CM,O_R,O_EC,Q,2,Wfwn_cloud);
      test_common::assert_near(Wfwn_cloud,Wexact,1e-2);
    }
  };
  // FWN clouds using barycenters won't work well for very coarse models like the cube
  test_common::run_test_cases<std::string>(
    {"bunny.off", "elephant.off", "hemisphere.obj"},
    test_case);
}
