#include <test_common.h>
#include <igl/linprog.h>

TEST_CASE("linprog: 2D-inequality", "[igl]" )
{
  Eigen::VectorXd f(2); f << -2, -1;
  Eigen::VectorXd b(3); b << 8, 14, 4;
  Eigen::MatrixXd A(3, 2); A << 2, -1, 1, 2, -1, 1;
  Eigen::MatrixXd B = Eigen::MatrixXd(0,2);
  Eigen::VectorXd c = Eigen::VectorXd(0);
  Eigen::VectorXd x;
  igl::linprog(f, A, b, B, c, x);
  Eigen::VectorXd x_correct(2); x_correct<<6,4;
  test_common::assert_near( x, x_correct, 1e-10);
}

TEST_CASE("linprog: 2D-inequality+2-equality", "[igl]" )
{
  Eigen::VectorXd f(2); f << -2, -1;
  Eigen::VectorXd b(3); b << 8, 14, 4;
  Eigen::MatrixXd A(3, 2); A << 2, -1, 1, 2, -1, 1;
  Eigen::MatrixXd B = Eigen::MatrixXd(2,2);B<<1,0,0,1;
  Eigen::VectorXd c = Eigen::VectorXd(2);c<<4,4;
  Eigen::VectorXd x;
  igl::linprog(f, A, b, B, c, x);
  Eigen::VectorXd x_correct(2); x_correct<<4,4;
  test_common::assert_near( x, x_correct, 1e-10);
}

TEST_CASE("linprog: 2D-inequality+1-equality", "[igl]" )
{
  Eigen::VectorXd f(2); f << -2, -1;
  Eigen::VectorXd b(3); b << 8, 14, 4;
  Eigen::MatrixXd A(3, 2); A << 2, -1, 1, 2, -1, 1;
  Eigen::MatrixXd B = Eigen::MatrixXd(1,2);B<<1,0;
  Eigen::VectorXd c = Eigen::VectorXd(1);c<<4;
  Eigen::VectorXd x;
  igl::linprog(f, A, b, B, c, x);
  Eigen::VectorXd x_correct(2); x_correct<<4,5;
  test_common::assert_near( x, x_correct, 1e-10);
}
