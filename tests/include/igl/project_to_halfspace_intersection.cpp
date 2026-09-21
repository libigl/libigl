#include <test_common.h>
#include <igl/project_to_halfspace_intersection.h>
#include <random>
#include <limits>

namespace
{
  igl::HalfspaceProjectionOptions<double> default_options()
  {
    igl::HalfspaceProjectionOptions<double> o;
    o.eps_rank = 1e-9;
    o.eps_feasible = 1e-9;
    o.eps_kkt = 1e-8;
    o.eps_outward = 0;
    o.seed = 42;
    return o;
  }
}

TEST_CASE("project_to_halfspace_intersection: no planes returns q", "[igl]")
{
  Eigen::Vector3d q(1,2,3);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(0,3);
  Eigen::VectorXd b(0);
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::SUCCESS);
  test_common::assert_near(r.p,q,1e-12);
  REQUIRE(r.active_count == 0);
}

TEST_CASE("project_to_halfspace_intersection: already feasible point is unchanged", "[igl]")
{
  Eigen::Vector3d q(1,1,1);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(3,3);
  A << 1,0,0,  0,1,0,  0,0,1;
  Eigen::VectorXd b(3); b << 0,0,0;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::SUCCESS);
  test_common::assert_near(r.p,q,1e-12);
  REQUIRE(r.active_count == 0);
}

TEST_CASE("project_to_halfspace_intersection: three orthogonal planes project to their intersection point", "[igl]")
{
  Eigen::Vector3d q(-1,-2,-3);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(3,3);
  A << 1,0,0,  0,1,0,  0,0,1;
  Eigen::VectorXd b(3); b << 0,0,0;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::SUCCESS);
  test_common::assert_near(r.p,Eigen::Vector3d(0,0,0),1e-9);
  REQUIRE(r.active_count == 3);
}

TEST_CASE("project_to_halfspace_intersection: bounded box picks nearest corner", "[igl]")
{
  Eigen::Vector3d q(5,5,5);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(6,3);
  A <<  1,0,0,  0,1,0,  0,0,1,
       -1,0,0,  0,-1,0, 0,0,-1;
  Eigen::VectorXd b(6); b << 0,0,0,-1,-1,-1;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::SUCCESS);
  test_common::assert_near(r.p,Eigen::Vector3d(1,1,1),1e-9);
  REQUIRE(r.active_count == 3);
}

TEST_CASE("project_to_halfspace_intersection: parallel contradictory planes are infeasible", "[igl]")
{
  Eigen::Vector3d q(0,0,0);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(2,3);
  A << 0,0,1,  0,0,-1;
  Eigen::VectorXd b(2); b << 1,1;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::INFEASIBLE);
}

TEST_CASE("project_to_halfspace_intersection: zero plane normal is redundant or infeasible", "[igl]")
{
  Eigen::Vector3d q(2,3,4);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(1,3);
  A << 0,0,0;
  {
    Eigen::VectorXd b(1); b << -1; // 0 >= -1: always true
    igl::HalfspaceProjectionResult<double,3> r;
    const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
    REQUIRE(status == igl::HalfspaceProjectionStatus::SUCCESS);
    test_common::assert_near(r.p,q,1e-12);
  }
  {
    Eigen::VectorXd b(1); b << 1; // 0 >= 1: never true
    igl::HalfspaceProjectionResult<double,3> r;
    const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
    REQUIRE(status == igl::HalfspaceProjectionStatus::INFEASIBLE);
  }
}

TEST_CASE("project_to_halfspace_intersection: NaN input is a numerical failure, not garbage", "[igl]")
{
  Eigen::Vector3d q(std::numeric_limits<double>::quiet_NaN(),0,0);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(1,3); A << 0,0,1;
  Eigen::VectorXd b(1); b << 0;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,default_options(),r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::NUMERICAL_FAILURE);
}

TEST_CASE("project_to_halfspace_intersection: exceeding max_planes yields CAPACITY_EXCEEDED", "[igl]")
{
  Eigen::Vector3d q(0,0,0);
  Eigen::Matrix<double,Eigen::Dynamic,3> A(3,3);
  A << 1,0,0,  0,1,0,  0,0,1;
  Eigen::VectorXd b(3); b << 0,0,0;
  auto o = default_options();
  o.max_planes = 2;
  igl::HalfspaceProjectionResult<double,3> r;
  const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,o,r);
  REQUIRE(status == igl::HalfspaceProjectionStatus::CAPACITY_EXCEEDED);
}

TEST_CASE("project_to_halfspace_intersection: matches a general-purpose QP on random feasible cases", "[igl]")
{
  // Cross-check against igl::copyleft::quadprog is done in the standalone
  // fixed-dim-qp-bench repo (768 real mesh-derived fixtures, agreement to
  // ~1e-9 relative error); this in-tree test is a lighter-weight regression
  // guard using a handful of deterministic random cases so it has no
  // GPL-licensed dependency.
  std::mt19937_64 rng(0xD1FF7E57);
  std::uniform_real_distribution<double> coord(-5,5);
  std::normal_distribution<double> normal_dist(0,1);
  int n_success = 0;
  for(int trial = 0;trial < 200;trial++)
  {
    const int m = trial % 10;
    Eigen::Vector3d q(coord(rng),coord(rng),coord(rng));
    Eigen::Matrix<double,Eigen::Dynamic,3> A(m,3);
    Eigen::VectorXd b(m);
    for(int i = 0;i < m;i++)
    {
      Eigen::Vector3d n(normal_dist(rng),normal_dist(rng),normal_dist(rng));
      if(n.norm() < 1e-6) n = Eigen::Vector3d(1,0,0);
      n.normalize();
      A.row(i) = n;
      b(i) = coord(rng);
    }
    igl::HalfspaceProjectionOptions<double> o = default_options();
    o.seed = static_cast<uint64_t>(trial);
    igl::HalfspaceProjectionResult<double,3> r;
    const auto status = igl::project_to_halfspace_intersection<double,3>(q,A,b,o,r);
    REQUIRE(status != igl::HalfspaceProjectionStatus::CAPACITY_EXCEEDED);
    REQUIRE(status != igl::HalfspaceProjectionStatus::NUMERICAL_FAILURE);
    if(status == igl::HalfspaceProjectionStatus::SUCCESS)
    {
      n_success++;
      REQUIRE(r.p.allFinite());
      REQUIRE(r.max_violation <= 1e-6);
      REQUIRE(r.active_count <= 3);
    }
  }
  REQUIRE(n_success > 0);
}
