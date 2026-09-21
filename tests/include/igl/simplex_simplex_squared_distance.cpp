#include <test_common.h>
#include <igl/simplex_simplex_squared_distance.h>
#include <random>
#include <vector>

namespace
{
  // Check the reported squared distance against `expected` and sanity check
  // that the returned barycentric coordinates are valid and reproduce it.
  void check(
    const Eigen::MatrixXd & V1,
    const Eigen::MatrixXd & V2,
    const double expected,
    const double margin = 1e-12)
  {
    double sqrdDist;
    Eigen::VectorXd B1,B2;
    igl::simplex_simplex_squared_distance(V1,V2,sqrdDist,B1,B2);
    REQUIRE (sqrdDist == Approx(expected).margin(margin));
    REQUIRE (B1.size() == V1.rows());
    REQUIRE (B2.size() == V2.rows());
    REQUIRE (B1.sum() == Approx(1.0).margin(1e-12));
    REQUIRE (B2.sum() == Approx(1.0).margin(1e-12));
    REQUIRE (B1.minCoeff() >= -1e-12);
    REQUIRE (B2.minCoeff() >= -1e-12);
    const Eigen::RowVectorXd c1 = B1.transpose()*V1;
    const Eigen::RowVectorXd c2 = B2.transpose()*V2;
    REQUIRE ((c1-c2).squaredNorm() == Approx(expected).margin(margin));
  }

  // Exhaustive reference implementation: over every pair of nonempty faces
  // (subsets of corners) solve the affine least-squares problem and keep the
  // minimum over those whose barycentric weights are all nonnegative. No
  // pruning, no recursion.
  double exhaustive(const Eigen::MatrixXd & A, const Eigen::MatrixXd & B)
  {
    const int na = A.rows();
    const int nb = B.rows();
    const int dim = A.cols();
    double best = std::numeric_limits<double>::infinity();
    for(int ma = 1;ma<(1<<na);ma++)
    {
      std::vector<int> I;
      for(int i = 0;i<na;i++) { if(ma&(1<<i)) { I.push_back(i); } }
      for(int mb = 1;mb<(1<<nb);mb++)
      {
        std::vector<int> J;
        for(int i = 0;i<nb;i++) { if(mb&(1<<i)) { J.push_back(i); } }
        const int m1 = I.size();
        const int m2 = J.size();
        Eigen::VectorXd w1 = Eigen::VectorXd::Zero(na);
        Eigen::VectorXd w2 = Eigen::VectorXd::Zero(nb);
        if(m1 == 1 && m2 == 1)
        {
          w1(I[0]) = 1;
          w2(J[0]) = 1;
        }else
        {
          Eigen::MatrixXd M(dim,(m1-1)+(m2-1));
          for(int i = 1;i<m1;i++)
          {
            M.col(i-1) = (A.row(I[i])-A.row(I[0])).transpose();
          }
          for(int i = 1;i<m2;i++)
          {
            M.col(m1-1+i-1) = -(B.row(J[i])-B.row(J[0])).transpose();
          }
          const Eigen::VectorXd rhs = (B.row(J[0])-A.row(I[0])).transpose();
          const Eigen::VectorXd x = M.completeOrthogonalDecomposition().solve(rhs);
          double s1 = 0;
          double s2 = 0;
          for(int i = 1;i<m1;i++) { w1(I[i]) = x(i-1); s1 += x(i-1); }
          w1(I[0]) = 1-s1;
          for(int i = 1;i<m2;i++) { w2(J[i]) = x(m1-1+i-1); s2 += x(m1-1+i-1); }
          w2(J[0]) = 1-s2;
        }
        if(w1.minCoeff() < -1e-12 || w2.minCoeff() < -1e-12) { continue; }
        best = std::min(best,
          ((w1.transpose()*A)-(w2.transpose()*B)).squaredNorm());
      }
    }
    return best;
  }
}

TEST_CASE("simplex_simplex_squared_distance: point-point", "[igl]")
{
  Eigen::MatrixXd V1(1,3); V1 << 0,0,0;
  Eigen::MatrixXd V2(1,3); V2 << 1,2,2;
  check(V1,V2,9.0);
}

TEST_CASE("simplex_simplex_squared_distance: point-segment", "[igl]")
{
  Eigen::MatrixXd V2(2,2); V2 << -5,0, 5,0;
  {
    // Perpendicular foot lands inside the segment.
    Eigen::MatrixXd V1(1,2); V1 << 0,1;
    check(V1,V2,1.0);
  }
  {
    // Perpendicular foot lands outside: closest point is an endpoint.
    Eigen::MatrixXd V1(1,2); V1 << 10,1;
    check(V1,V2,26.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: segment-segment", "[igl]")
{
  {
    // Parallel.
    Eigen::MatrixXd V1(2,3); V1 << 0,0,0, 1,0,0;
    Eigen::MatrixXd V2(2,3); V2 << 0,1,0, 1,1,0;
    check(V1,V2,1.0);
  }
  {
    // Skew.
    Eigen::MatrixXd V1(2,3); V1 << -1,0,0, 1,0,0;
    Eigen::MatrixXd V2(2,3); V2 << 0,-1,1, 0,1,1;
    check(V1,V2,1.0);
  }
  {
    // Collinear and disjoint.
    Eigen::MatrixXd V1(2,3); V1 << 0,0,0, 1,0,0;
    Eigen::MatrixXd V2(2,3); V2 << 3,0,0, 5,0,0;
    check(V1,V2,4.0);
  }
  {
    // Collinear and overlapping.
    Eigen::MatrixXd V1(2,3); V1 << 0,0,0, 2,0,0;
    Eigen::MatrixXd V2(2,3); V2 << 1,0,0, 5,0,0;
    check(V1,V2,0.0);
  }
  {
    // Crossing in 2D.
    Eigen::MatrixXd V1(2,2); V1 << -1,0, 1,0;
    Eigen::MatrixXd V2(2,2); V2 << 0,-1, 0,1;
    check(V1,V2,0.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: point-triangle", "[igl]")
{
  Eigen::MatrixXd V2(3,3); V2 << 0,0,0, 1,0,0, 0,1,0;
  {
    // Projects into the interior.
    Eigen::MatrixXd V1(1,3); V1 << 0.25,0.25,3;
    check(V1,V2,9.0);
  }
  {
    // Projects outside: closest point is a corner.
    Eigen::MatrixXd V1(1,3); V1 << -1,-1,0;
    check(V1,V2,2.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: triangle-triangle", "[igl]")
{
  {
    // Parallel planes.
    Eigen::MatrixXd V1(3,3); V1 << 0,0,0, 1,0,0, 0,1,0;
    Eigen::MatrixXd V2(3,3); V2 << 0,0,2, 1,0,2, 0,1,2;
    check(V1,V2,4.0);
  }
  {
    // Intersecting.
    Eigen::MatrixXd V1(3,3); V1 << -1,-1,0, 3,-1,0, -1,3,0;
    Eigen::MatrixXd V2(3,3); V2 << 0.2,0.2,-1, 0.2,0.2,1, 1,1,1;
    check(V1,V2,0.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: tet", "[igl]")
{
  Eigen::MatrixXd V1(4,3); V1 << 0,0,0, 1,0,0, 0,1,0, 0,0,1;
  {
    Eigen::MatrixXd V2(4,3); V2 << 0,0,5, 1,0,5, 0,1,5, 0,0,6;
    check(V1,V2,16.0);
  }
  {
    // Point strictly inside.
    Eigen::MatrixXd V2(1,3); V2 << 0.2,0.2,0.2;
    check(V1,V2,0.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: degenerate", "[igl]")
{
  {
    // "Triangle" with a repeated corner is really a segment.
    Eigen::MatrixXd V1(3,3); V1 << 0,0,0, 1,0,0, 1,0,0;
    Eigen::MatrixXd V2(2,3); V2 << 0.5,2,0, 0.5,3,0;
    check(V1,V2,4.0);
  }
  {
    // All corners identical.
    Eigen::MatrixXd V1(3,3); V1 << 1,1,1, 1,1,1, 1,1,1;
    Eigen::MatrixXd V2(1,3); V2 << 1,1,4;
    check(V1,V2,9.0);
  }
}

TEST_CASE("simplex_simplex_squared_distance: output-types", "[igl]")
{
  {
    // Row-vector outputs.
    Eigen::MatrixXd V1(2,3); V1 << 0,0,0, 1,0,0;
    Eigen::MatrixXd V2(1,3); V2 << 0.5,1,0;
    double sqrdDist;
    Eigen::RowVectorXd B1,B2;
    igl::simplex_simplex_squared_distance(V1,V2,sqrdDist,B1,B2);
    REQUIRE (sqrdDist == Approx(1.0).margin(1e-12));
    REQUIRE (B1.rows() == 1);
    REQUIRE (B1.cols() == 2);
    REQUIRE (B1(0) == Approx(0.5).margin(1e-12));
    REQUIRE (B1(1) == Approx(0.5).margin(1e-12));
  }
  {
    // Fixed-size float inputs and outputs.
    Eigen::Matrix<float,2,3> V1; V1 << 0,0,0, 1,0,0;
    Eigen::Matrix<float,3,3> V2; V2 << 0,2,0, 1,2,0, 0,3,0;
    float sqrdDist;
    Eigen::Matrix<float,2,1> B1;
    Eigen::Matrix<float,3,1> B2;
    igl::simplex_simplex_squared_distance(V1,V2,sqrdDist,B1,B2);
    REQUIRE (sqrdDist == Approx(4.0f).margin(1e-6));
  }
  {
    // Mixed scalar types.
    Eigen::MatrixXd V1(1,3); V1 << 0,0,0;
    Eigen::MatrixXf V2(1,3); V2 << 0,0,3;
    double sqrdDist;
    Eigen::VectorXd B1,B2;
    igl::simplex_simplex_squared_distance(V1,V2,sqrdDist,B1,B2);
    REQUIRE (sqrdDist == Approx(9.0).margin(1e-12));
  }
}

TEST_CASE("simplex_simplex_squared_distance: fixed-size matches dynamic", "[igl]")
{
  // Fixed-size inputs take a separate, compile-time-unrolled code path. It must
  // agree with the dynamically sized one.
  std::mt19937 gen(1337);
  std::uniform_real_distribution<double> U(-1,1);
  for(int trial = 0;trial<500;trial++)
  {
    // Half the time push the simplices apart, half the time let them overlap.
    const double shift = (gen()%2) ? 1.5 : 0.0;
    Eigen::Matrix3d T1,T2;
    Eigen::Matrix<double,4,3> E1,E2;
    for(int i = 0;i<3;i++)
    {
      for(int j = 0;j<3;j++)
      {
        T1(i,j) = U(gen);
        T2(i,j) = U(gen) + (j==0?shift:0.0);
      }
    }
    for(int i = 0;i<4;i++)
    {
      for(int j = 0;j<3;j++)
      {
        E1(i,j) = U(gen);
        E2(i,j) = U(gen) + (j==0?shift:0.0);
      }
    }
    // Occasionally degenerate.
    if(trial%5 == 0) { T1.row(2) = T1.row(0); }
    if(trial%7 == 0) { E2.row(3) = 0.5*(E2.row(0)+E2.row(1)); }
    {
      double fixed_sqrd,dynamic_sqrd;
      Eigen::Vector3d fixed_B1,fixed_B2;
      Eigen::VectorXd dynamic_B1,dynamic_B2;
      igl::simplex_simplex_squared_distance(
        T1,T2,fixed_sqrd,fixed_B1,fixed_B2);
      const Eigen::MatrixXd D1 = T1, D2 = T2;
      igl::simplex_simplex_squared_distance(
        D1,D2,dynamic_sqrd,dynamic_B1,dynamic_B2);
      REQUIRE (fixed_sqrd == Approx(dynamic_sqrd).margin(1e-9));
      REQUIRE (fixed_B1.sum() == Approx(1.0).margin(1e-12));
      REQUIRE (fixed_B1.minCoeff() >= -1e-12);
      REQUIRE (((fixed_B1.transpose()*T1)-(fixed_B2.transpose()*T2)).squaredNorm()
        == Approx(fixed_sqrd).margin(1e-9));
    }
    {
      double fixed_sqrd,dynamic_sqrd;
      Eigen::Vector4d fixed_B1,fixed_B2;
      Eigen::VectorXd dynamic_B1,dynamic_B2;
      igl::simplex_simplex_squared_distance(
        E1,E2,fixed_sqrd,fixed_B1,fixed_B2);
      const Eigen::MatrixXd D1 = E1, D2 = E2;
      igl::simplex_simplex_squared_distance(
        D1,D2,dynamic_sqrd,dynamic_B1,dynamic_B2);
      REQUIRE (fixed_sqrd == Approx(dynamic_sqrd).margin(1e-9));
      REQUIRE (fixed_B2.sum() == Approx(1.0).margin(1e-12));
      REQUIRE (fixed_B2.minCoeff() >= -1e-12);
    }
  }
}

TEST_CASE("simplex_simplex_squared_distance: exhaustive", "[igl]")
{
  // Compare the pruned recursion against an unpruned enumeration of every
  // face-pair subproblem, over random (and deliberately degenerate) simplices.
  std::mt19937 gen(42);
  std::uniform_real_distribution<double> U(-1,1);
  for(int trial = 0;trial<500;trial++)
  {
    const int dim = 1 + trial%3;
    const int na = 1 + (int)(gen()%4);
    const int nb = 1 + (int)(gen()%4);
    Eigen::MatrixXd V1(na,dim),V2(nb,dim);
    // Half the time push the simplices apart, half the time let them overlap.
    const double shift = (gen()%2) ? 2.5 : 0.0;
    for(int i = 0;i<na;i++) { for(int j = 0;j<dim;j++) { V1(i,j) = U(gen); } }
    for(int i = 0;i<nb;i++)
    {
      for(int j = 0;j<dim;j++) { V2(i,j) = U(gen) + (j==0?shift:0.0); }
    }
    // Sprinkle in degeneracies: duplicate and collinear corners.
    if(trial%5 == 0 && na>1) { V1.row(na-1) = V1.row(0); }
    if(trial%7 == 0 && nb>2) { V2.row(nb-1) = 0.5*(V2.row(0)+V2.row(1)); }
    // And exact coincidences by snapping to a coarse grid.
    if(trial%11 == 0)
    {
      V1 = (V1.array()*4).round()/4;
      V2 = (V2.array()*4).round()/4;
    }
    check(V1,V2,exhaustive(V1,V2),1e-9);
  }
}
