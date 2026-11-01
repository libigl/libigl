#include <test_common.h>
#include <igl/active_set.h>

TEST_CASE("active_set: empty initial guess", "[igl]")
{
  Eigen::SparseMatrix<double> A(3,3);
  A.setIdentity();
  const Eigen::VectorXd B = Eigen::VectorXd::Zero(3);
  const Eigen::VectorXd lx = Eigen::VectorXd::Constant(3,-10);
  const Eigen::VectorXd ux = Eigen::VectorXd::Constant(3,10);
  for(const bool provide_initial_guess : {false,true})
  {
    CAPTURE(provide_initial_guess);
    Eigen::VectorXd Z;
    if(provide_initial_guess)
    {
      Z = B;
    }
    const auto status = igl::active_set(
      A,B,Eigen::VectorXi(),Eigen::VectorXd(),
      Eigen::SparseMatrix<double>(),Eigen::VectorXd(),
      Eigen::SparseMatrix<double>(),Eigen::VectorXd(),
      lx,ux,igl::active_set_params(),Z);
    REQUIRE(status == igl::SOLVER_STATUS_CONVERGED);
    test_common::assert_eq(Z,B);
  }
}

TEST_CASE("active_set: nonzero solution", "[igl]")
{
  Eigen::SparseMatrix<double> A(3,3);
  A.insert(0,0) = 2;
  A.insert(1,1) = 3;
  A.insert(2,2) = 4;
  Eigen::VectorXd B(3);
  B << -2,6,-12;
  Eigen::VectorXi known;
  Eigen::VectorXd Y,Beq,Bieq;
  Eigen::SparseMatrix<double> Aeq(0,3),Aieq(0,3);
  Eigen::VectorXd lx = Eigen::VectorXd::Constant(3,-10);
  Eigen::VectorXd ux = Eigen::VectorXd::Constant(3,10);
  Eigen::VectorXd expected(3);
  expected << 1,-2,3;
  igl::active_set_params params;
  params.Auu_pd = true;

  for(const bool provide_initial_guess : {false,true})
  {
    CAPTURE(provide_initial_guess);
    Eigen::VectorXd Z;
    if(provide_initial_guess)
    {
      Z = Eigen::VectorXd::Constant(3,0.5);
    }
    const auto status = igl::active_set(
      A,B,known,Y,Aeq,Beq,Aieq,Bieq,lx,ux,params,Z);
    REQUIRE(status == igl::SOLVER_STATUS_CONVERGED);
    test_common::assert_near(Z,expected,1e-12);
  }
}

TEST_CASE("active_set: active constraints", "[igl]")
{
  Eigen::SparseMatrix<double> A(3,3);
  A.setIdentity();
  Eigen::VectorXd B(3);
  B << 2,-3,-4;
  Eigen::VectorXi known;
  Eigen::VectorXd Y,Beq;
  Eigen::SparseMatrix<double> Aeq(0,3),Aieq(1,3);
  Aieq.insert(0,2) = 1;
  Eigen::VectorXd Bieq = Eigen::VectorXd::Constant(1,2);
  Eigen::VectorXd lx(3),ux(3),expected(3);
  lx << 0,-10,-10;
  ux << 10,1,10;
  expected << 0,1,2;
  igl::active_set_params params;
  params.Auu_pd = true;

  for(const bool provide_initial_guess : {false,true})
  {
    CAPTURE(provide_initial_guess);
    Eigen::VectorXd Z;
    if(provide_initial_guess)
    {
      Z = Eigen::VectorXd::Zero(3);
    }
    const auto status = igl::active_set(
      A,B,known,Y,Aeq,Beq,Aieq,Bieq,lx,ux,params,Z);
    REQUIRE(status == igl::SOLVER_STATUS_CONVERGED);
    test_common::assert_near(Z,expected,1e-12);
    REQUIRE((Z.array() >= lx.array()-1e-12).all());
    REQUIRE((Z.array() <= ux.array()+1e-12).all());
    REQUIRE(((Aieq*Z-Bieq).array() <= 1e-12).all());
  }
}

TEST_CASE("active_set: activate constraints after a small step", "[igl]")
{
  const double delta = 1e-8;
  Eigen::SparseMatrix<double> A(1,1);
  A.setIdentity();
  Eigen::VectorXd B = Eigen::VectorXd::Constant(1,-1-delta);
  Eigen::VectorXi known;
  Eigen::VectorXd Y,Beq,Bieq;
  Eigen::SparseMatrix<double> Aeq(0,1),Aieq(0,1);
  Eigen::VectorXd lx = Eigen::VectorXd::Constant(1,-10);
  Eigen::VectorXd ux = Eigen::VectorXd::Constant(1,10);
  Eigen::VectorXd Z = Eigen::VectorXd::Constant(1,1-delta);
  igl::active_set_params params;
  params.Auu_pd = true;

  SECTION("lower bound")
  {
    lx(0) = 1;
    B(0) = -1+delta;
    Z(0) = 1+delta;
  }
  SECTION("upper bound")
  {
    ux(0) = 1;
  }
  SECTION("linear inequality")
  {
    Aieq.resize(1,1);
    Aieq.insert(0,0) = 1;
    Bieq = Eigen::VectorXd::Constant(1,1);
  }

  // The unconstrained step is small enough for convergence, but violates
  // a constraint that must be enforced by another solve
  REQUIRE(4*delta*delta < params.solution_diff_threshold);
  const auto status = igl::active_set(
    A,B,known,Y,Aeq,Beq,Aieq,Bieq,lx,ux,params,Z);
  REQUIRE(status == igl::SOLVER_STATUS_CONVERGED);
  REQUIRE(std::abs(Z(0)-1) < 1e-12);
}

TEST_CASE("active_set: release constraints after a small step", "[igl]")
{
  const double delta = 1e-8;
  Eigen::SparseMatrix<double> A(1,1);
  A.setIdentity();
  Eigen::VectorXd B = Eigen::VectorXd::Constant(1,-0.5);
  Eigen::VectorXi known;
  Eigen::VectorXd Y,Beq,Bieq;
  Eigen::SparseMatrix<double> Aeq(0,1),Aieq(0,1);
  Eigen::VectorXd lx = Eigen::VectorXd::Constant(1,-10);
  Eigen::VectorXd ux = Eigen::VectorXd::Constant(1,10);
  Eigen::VectorXd Z = Eigen::VectorXd::Constant(1,1+delta);
  igl::active_set_params params;
  params.Auu_pd = true;

  SECTION("lower bound")
  {
    lx(0) = 1;
    B(0) = -2;
    Z(0) = 1-delta;
  }
  SECTION("upper bound")
  {
    ux(0) = 1;
  }
  SECTION("linear inequality")
  {
    Aieq.resize(1,1);
    Aieq.insert(0,0) = 1;
    Bieq = Eigen::VectorXd::Constant(1,1);
  }

  // The first solve barely changes the infeasible initial guess
  // Its negative multiplier releases the constraint, so another solve is needed
  REQUIRE(delta*delta < params.solution_diff_threshold);
  const auto status = igl::active_set(
    A,B,known,Y,Aeq,Beq,Aieq,Bieq,lx,ux,params,Z);
  REQUIRE(status == igl::SOLVER_STATUS_CONVERGED);
  REQUIRE(std::abs(Z(0)+B(0)) < 1e-12);
}
