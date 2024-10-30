#include <test_common.h>
#include <igl/triangle/refine.h>
#include <igl/triangle/triangulate.h>
#include <igl/matlab_format.h>
#include <igl/find.h>
#include <igl/barycenter.h>
#include <igl/winding_number.h>
#include <igl/doublearea.h>
#include <iostream>

TEST_CASE("refine", "[igl][triangle]") {
  Eigen::MatrixXd V(3,2);
  V << 
    0,0,
    1,0,
    0,1;
  Eigen::MatrixXi F(1,3);
  F<<0,1,2;
  std::string flags = "a0.0625Q";
  Eigen::MatrixXd V2;
  Eigen::MatrixXi F2;
  Eigen::MatrixXi E;
  igl::triangle::refine(V,E,F,flags,V2,F2);
  REQUIRE(F2.rows()>1);
  Eigen::VectorXd A2;
  igl::doublearea(V2,F2,A2);
  REQUIRE((A2.array()<2*0.0625).all());
}


TEST_CASE("refine-inner-pinch", "[igl][triangle]") {
  // Mesh an outline that looks like:
  //
  // o-----------------o
  // |                 |
  // |   o         o   |
  // |   |\       /|   |
  // |   | \     / |   |
  // |   |  \   /  |   |
  // |   |   \ /   |   |
  // |   |    o    |   |
  // |   o---------o   |
  // |                 |
  // o-----------------o
  //
  // If we mesh the whole convex hull and then remove the inner cavity we'll get
  // over-refinement near the sharp corner.
  //
  // We could tell triangle a point inside this cavity but that's not often easy
  // to do robustly/for broken input. Instead we can use a trivial mesh of the
  // whole domain, remove the cavity with the winding number and then refine
  // what's left.
  //
  Eigen::MatrixXd V(9,2);
  V << 
    0,0,
    3,0,
    3,3,
    0,3,
    1,1,
    2,1,
    2,2,
    1.5,1.00000001,
    1,2;
  Eigen::MatrixXi E(9,2);
  E <<
    0,1,
    1,2,
    2,3,
    3,0,
    8,7,
    7,6,
    6,5,
    5,4,
    4,8;

  Eigen::MatrixXd H(0,2);
  const auto filter  = [&V,&E](
    const Eigen::MatrixXd & Vc,
    Eigen::MatrixXi & Fc)
  {
    Eigen::MatrixXd BC;
    igl::barycenter(Vc,Fc,BC);
    Eigen::VectorXd W;
    // This is sadly O(|E| ⋅ |BC|) for 2D inputs.
    igl::winding_number(V,E,BC,W);
    auto keep = igl::find((W.array().abs() > 0.5).eval());
    Fc = Fc(keep,Eigen::all).eval();
  };

  Eigen::MatrixXd Vq;
  Eigen::MatrixXi Fq;
  igl::triangle::triangulate(V,E,H,"cq33a0.0625Q",Vq,Fq);
  filter(Vq,Fq);
  //std::cout<<igl::matlab_format(Vq,"Vq")<<std::endl;
  //std::cout<<igl::matlab_format_index(Fq,"Fq")<<std::endl;

  Eigen::MatrixXd Vc;
  Eigen::MatrixXi Fc;
  igl::triangle::triangulate(V,E,H,"cQ",Vc,Fc);
  REQUIRE(V.rows() == Vc.rows());
  REQUIRE(Fc.rows() == 12);
  filter(Vc,Fc);
  REQUIRE(Fc.rows() == 9);
  //std::cout<<igl::matlab_format(Vc,"Vc")<<std::endl;
  //std::cout<<igl::matlab_format_index(Fc,"Fc")<<std::endl;

  Eigen::MatrixXd Vr;
  Eigen::MatrixXi Fr;

  Eigen::MatrixXi _;
  igl::triangle::refine(Vc,_,Fc,"q33a0.0625Q",Vr,Fr);
  //std::cout<<igl::matlab_format(Vr,"Vr")<<std::endl;
  //std::cout<<igl::matlab_format_index(Fr,"Fr")<<std::endl;
  
  REQUIRE(Fr.rows() > Fc.rows());
  REQUIRE(Fr.rows() < 0.15 * Fq.rows());
}
