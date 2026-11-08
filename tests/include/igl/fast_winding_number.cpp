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
#include <Eigen/Geometry>
#include <cmath>
#include <functional>

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

namespace
{
// Varying face normals make first moments about each node's centroid nonzero.
void moment_soup(const int side, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
  const int count = side*side*side;
  V.resize(3*count,3);
  F.resize(count,3);
  for(int i = 0; i < count; ++i)
  {
    const Eigen::RowVector3d p(i%side,(i/side)%side,i/(side*side));
    V.row(3*i) = p;
    V.row(3*i+1) = p + Eigen::RowVector3d(.75,.125*(i%3),.125*(i%5));
    V.row(3*i+2) = p + Eigen::RowVector3d(.125*(i%7),.625,.125*(i%4));
    F.row(i) << 3*i,3*i+1,3*i+2;
  }
}
}

TEST_CASE("fast_winding_number: translated_second_moments", "[igl]")
{
  using namespace igl::FastWindingNumber::HDK_Sample;
  for(const int side : {4,8})
  {
    CAPTURE(side);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    moment_soup(side,V,F);
    std::vector<UT_Vector3T<float>> positions(V.rows());
    std::vector<int> triangles(F.size());
    std::vector<UT::Box<float,3>> boxes(F.rows());
    for(int i = 0; i < V.rows(); ++i)
      for(int j = 0; j < 3; ++j) positions[i][j] = V(i,j);
    for(int i = 0; i < F.rows(); ++i)
    {
      boxes[i].initBounds(positions[F(i,0)]);
      for(int j = 0; j < 3; ++j)
      {
        triangles[3*i+j] = F(i,j);
        boxes[i].enlargeBounds(positions[F(i,j)]);
      }
    }
    // Rebuild the same BVH to identify the groups approximated at the root.
    UT_BVH<4> tree;
    tree.init<UT::BVH_Heuristic::BOX_AREA,float,3>(boxes.data(),F.rows());
    using Node = UT_BVH<4>::Node;
    const auto *nodes = tree.getNodes();
    std::vector<std::vector<int>> groups(4);
    std::function<void(unsigned,std::vector<int>&)> collect;
    collect = [&](unsigned child, std::vector<int> &items)
    {
      if(child == Node::EMPTY) return;
      if(!Node::isInternal(child)) { items.push_back(child); return; }
      for(const auto next : nodes[Node::getInternalNum(child)].child)
        collect(next,items);
    };
    for(int i = 0; i < 4; ++i)
    {
      collect(nodes[0].child[i],groups[i]);
      // A triangle's first moments vanish; deeper groups exercise their shift.
      REQUIRE(groups[i].size() > 4);
    }
    UT_SolidAngle<float,float> order2(F.rows(),triangles.data(),V.rows(),positions.data(),2);
    for(const Eigen::Vector3d direction : {Eigen::Vector3d(6,-2,7),
      Eigen::Vector3d(-4,7,5),Eigen::Vector3d(7,8,-3)})
    {
      const Eigen::Vector3d query = direction*(side/4.);
      double expected = 0;
      for(const auto &group : groups)
      {
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        double area = 0;
        for(const int i : group)
        {
          const Eigen::Vector3d a=V.row(F(i,0)), b=V.row(F(i,1)), c=V.row(F(i,2));
          const double weight = (b-a).cross(c-a).norm();
          center += weight*(a+b+c)/3.;
          area += weight;
        }
        center /= area;
        const Eigen::Vector3d q = query-center;
        const double r2 = q.squaredNorm(), r5 = r2*r2*std::sqrt(r2);
        for(const int i : group)
        {
          const Eigen::Vector3d a=V.row(F(i,0)), b=V.row(F(i,1)), c=V.row(F(i,2));
          const Eigen::Vector3d N = .5*(b-a).cross(c-a);
          // Three edge midpoints integrate the second-order Taylor polynomial
          // exactly, without using the moments or their translation formula.
          for(int j = 0; j < 3; ++j)
          {
            const Eigen::Vector3d u = .5*(V.row(F(i,j))+V.row(F(i,(j+1)%3))).transpose()-center;
            const double qu = q.dot(u), nq = N.dot(q);
            const double nu = N.dot(u);
            expected += ((nu-nq)*r2-3*nq*qu
              +3*nu*qu+1.5*nq*u.squaredNorm()-7.5*nq*qu*qu/r2)/(3*r5);
          }
        }
      }
      UT_Vector3T<float> p;
      for(int j = 0; j < 3; ++j) p[j] = query[j];
      // Zero accuracy_scale evaluates the root's children without descending.
      // Comparing the same Taylor polynomial avoids any truncation remainder.
      const double actual = order2.computeSolidAngle(p,0);
      CAPTURE(query,expected,actual);
      CHECK(std::abs(actual-expected) < 1e-7);
    }
  }
}

TEST_CASE("fast_winding_number: translated_soup_moments", "[igl]")
{
  // Vary both positions and normals so child first moments are nonzero.
  // The 64 triangles require several levels in the four-way hierarchy.
  Eigen::MatrixXd V(3*64,3);
  Eigen::MatrixXi F(64,3);
  for(int i = 0; i < 64; ++i)
  {
    const double x = (i%4)/3.0-0.2;
    const double y = ((i/4)%4)/3.0-0.3;
    const double z = (i/16)/3.0-0.1;
    const Eigen::RowVector3d center(
      x+0.3*y-0.2*z, y-0.2*x+0.3*z, z+0.2*x-0.1*y);
    const Eigen::RowVector3d u(0.2,0,-y*y);
    const Eigen::RowVector3d v(0,0.2,-z*z);
    V.row(3*i) = center-(u+v)/3.0;
    V.row(3*i+1) = V.row(3*i)+u;
    V.row(3*i+2) = V.row(3*i)+v;
    F.row(i) << 3*i,3*i+1,3*i+2;
  }
  // The soup implementation stores float vertices. Use the same geometry
  // for the double-precision direct reference.
  V = V.cast<float>().cast<double>().eval();
  igl::FastWindingNumberBVH bvh;
  igl::fast_winding_number(V,F,2,bvh);
  Eigen::MatrixXd Q(3,3);
  Q << 8,16,24, -16,8,24, 16,-24,8;
  Eigen::VectorXd actual, expected;
  igl::fast_winding_number(bvh,2.0f,Q,actual);
  igl::winding_number(V,F,Q,expected);
  // Far-field queries use internal-node expansions. This tolerance allows
  // truncation and float rounding, but detects the missing linear terms.
  test_common::assert_near(actual,expected,2e-10);
}
