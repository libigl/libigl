#include <test_common.h>
#include <igl/AABB.h>
#include <igl/EPS.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/colon.h>
#include <igl/get_seconds.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/randperm.h>
#include <igl/read_triangle_mesh.h>
#include <iostream>

TEST_CASE("AABB: find_2d", "[igl]")
{
  Eigen::MatrixXd V(6,2);
  V << 
  0,0,
  1,0,
  0,1,
  2,1,
  2,2,
  1,2;
  Eigen::MatrixXi F(4,3);
  F<< 
  2,0,1,
  2,1,5,
  5,3,4,
  5,1,3;
  igl::AABB<Eigen::MatrixXd,2> tree;
  tree.init(V,F);
  Eigen::RowVector2d q(0.5,0.5);
  std::vector<int> r = tree.find(V,F,q);
  REQUIRE(r.size() == 2);
  REQUIRE(r[0] == 0);
  REQUIRE(r[1] == 1);
}

TEST_CASE("AABB: find_3d", "[igl]")
{
  Eigen::MatrixXd V(7,3);
  V << 0,0,1,
    1,0,1,
    0,1,1,
    2,1,1,
    2,2,1,
    1,2,1,
    0,0,0;

  Eigen::MatrixXi F(4,4);
  F << 
    0,1,2,6,
    1,3,2,6,
    3,4,5,6,
    3,5,1,6;

  igl::AABB<Eigen::MatrixXd,3> tree;
  tree.init(V,F);
  Eigen::RowVector3d q(0.5,0.5,1.0);
  std::vector<int> r = tree.find(V,F,q);
  REQUIRE(r.size() == 2);
  REQUIRE(r[0] == 0);
  REQUIRE(r[1] == 1);
}

TEST_CASE("AABB: insert", "[igl]")
{
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  for(int i = 0;i<3;i++)
  {
    igl::AABB<Eigen::MatrixXd, 3> * leaf = new igl::AABB<Eigen::MatrixXd, 3>();
    leaf->m_box = Eigen::AlignedBox<double, 3>(
      Eigen::RowVector3d(-i,-i,-i),
      Eigen::RowVector3d(i,i,i));
    auto * ret = tree->insert(leaf);
    tree = ret->root();
    REQUIRE(leaf->is_leaf());
  }
  REQUIRE(tree->is_root());
  delete tree;
}

TEST_CASE("AABB: dynamic", "[igl]")
{
  const int MAX_RUNS = 10;
  const auto test_case = [&MAX_RUNS](const std::string &param)
  {
    IGL_TICTOC_LAMBDA;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Load example mesh: GetParam() will be name of mesh file
    igl::read_triangle_mesh(test_common::data_path(param), V, F);
    // Make into soup
    V = V(Eigen::Map<Eigen::VectorXi>(F.data(),F.size()), Eigen::all).eval();
    F = Eigen::Map<Eigen::MatrixXi>(igl::colon<int>(0,V.rows()-1).data(),V.rows()/3,3).eval();
    Eigen::MatrixXd BC;
    igl::barycenter(V,F,BC);
    Eigen::VectorXi I_gt = igl::colon<int>(0,F.rows()-1);
    Eigen::VectorXd sqrD_gt = Eigen::VectorXd::Zero(F.rows());
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
    tree->init(V,F);
    tictoc();
    tree->squared_distance(V,F,BC,sqrD,I,C);
    const double t0 = tictoc();
    test_common::assert_eq(I,I_gt);
    test_common::assert_near(sqrD,sqrD_gt,1e-15);
    REQUIRE(tree->size() == 2*F.rows()-1);

    const double pad = igl::EPS<double>();
    const double h = igl::avg_edge_length(V,F);
    const auto leaves = tree->gather_leaves(F.rows());
    for(auto * leaf : leaves)
    {
      REQUIRE(leaf);
      REQUIRE(leaf->is_leaf());
    }
    // This is slow. Only test on small meshes
    if(F.rows() < 4000)
    {
      // Gather list of pointers to leaves and pad by ε to force rebuild
      tree = tree->pad(leaves,pad,2);
      REQUIRE(tree->is_root());
      REQUIRE(tree == tree->root());
      for(auto * leaf : leaves)
      {
        REQUIRE(leaf);
        REQUIRE(leaf->is_leaf());
      }
      tictoc();
      tree->squared_distance(V,F,BC,sqrD,I,C);
      const double t1 = tictoc();
      REQUIRE(t1 < 10*t0);
      test_common::assert_eq(I,I_gt);
      test_common::assert_near(sqrD,sqrD_gt,1e-15);
      REQUIRE(tree->size() == 2*F.rows()-1);
    }
#ifndef NDEBUG
    if(F.rows()>10000)
    {
      std::cerr<<"#ifndef NDEBUG: Skipping timing test."<<std::endl;
      delete tree;
      return;
    }
#endif


    Eigen::VectorXi RV;
    Eigen::VectorXi RF;
    // Perturb a small subset of the triangles
    {
      {
        igl::randperm(F.rows(),RF);
        RF = RF.topRows(std::min(12,(int)F.rows())).eval();
        RV.resize(RF.size()*3);
        RV << RF, RF.array()+F.rows(), RF.array()+2*F.rows();
      }
      Eigen::MatrixXd TF = 0.1*h*Eigen::MatrixXd::Random(RF.size(),3);
      Eigen::MatrixXd TV(RV.rows(),3);
      TV<<TF,TF,TF;
      V(RV,Eigen::all) += TV;
      igl::barycenter(V,F,BC);
    }
    const int qi = RF(0);

    double t_dynamic_tree;
    {
      tictoc();
      // update those perturbed triangles
      for(int i = 0;i<RF.size();i++)
      {
        tree = leaves[RF(i)]->update_primitive(V,F,pad)->root();
      }
      const double t_refit = tictoc();
      for(int r = 0;r<MAX_RUNS;r++)
      {
        tree->squared_distance(V,F,Eigen::MatrixXd(BC.row(qi)),sqrD,I,C);
      }
      const double t_query = tictoc()/MAX_RUNS;
      REQUIRE(I(0) == qi);
      t_dynamic_tree = t_refit+t_query;
    }
    double t_static_tree;
    {
      tictoc();
      for(int r = 0;r<MAX_RUNS;r++)
      {
        igl::point_mesh_squared_distance(
          Eigen::MatrixXd(BC.row(qi)),V,F,sqrD,I,C);
      }
      t_static_tree = tictoc()/MAX_RUNS;
      REQUIRE(I(0) == qi);
    }
    double t_brute_force;
    {
      tictoc();
      double min_dist;
      int min_i;
      for(int r = 0;r<MAX_RUNS;r++)
      {
        min_dist = std::numeric_limits<double>::infinity();
        min_i = -1;
        Eigen::RowVector3d q = Eigen::RowVector3d(BC.row(qi));
        for(int i = 0;i<F.rows();i++)
        {
          Eigen::RowVector3d c;
          double d;
          igl::point_simplex_squared_distance<3>(q,V,F,i,d,c);
          if(d < min_dist)
          {
            min_dist = d;
            min_i = i;
          }
        }
      }
      t_brute_force = tictoc()/MAX_RUNS;
      REQUIRE(min_i == qi);
    }

    // Only compare speeds on large meshes
    if(F.rows() > 300)
    { 
      REQUIRE(t_dynamic_tree < t_static_tree);
      REQUIRE(t_dynamic_tree < t_brute_force);
    }
    delete tree;
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}
