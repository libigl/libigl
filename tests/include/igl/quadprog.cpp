#include <test_common.h>
#include <igl/quadprog.h>
#include <igl/EPS.h>


TEST_CASE("quadprog: box3", "[igl]" )
{
  {
    Eigen::Matrix3d H = (Eigen::Matrix3d(3,3)<<
        0.240548455386281,   0.237314308102107,  0.0436993831501944,
        0.237314308102107,   0.326254049041854,-0.00021896091952631,
        0.0436993831501944,-0.00021896091952631,   0.171175681280756
        ).finished();
    Eigen::Vector3d f = (Eigen::Vector3d(3,1)<<
        0.222182612270718,
        0.503254616893693,
        -0.184619987497072
        ).finished();
    Eigen::Vector3d lb = (Eigen::Vector3d(3,1)<<
        -1,
        -1,
        -1
        ).finished();
    Eigen::Vector3d ub = (Eigen::Vector3d(3,1)<<
        1,
        1,
        1
        ).finished(); 
	// Windows needs template args spelled out
    Eigen::Vector3d x = igl::quadprog<double,3>(H,f,lb,ub);
    //std::cout<<igl::matlab_format(x,"x")<<std::endl;
    REQUIRE(abs(x(0)- -0.118760635036839)<1e-7);
    REQUIRE(abs(x(1)- -1)<1e-7);
    REQUIRE(abs(x(2)- +1)<1e-7);
  }
  {
    Eigen::Matrix3d H;
    H<<1,0,0,0,1,0,0,0,1;
    Eigen::Vector3d f( 0.5,-0.5,-0.5);
    Eigen::Vector3d lb(0,0,0);
    Eigen::Vector3d ub(1,1,1);
	// Windows needs template args spelled out
    Eigen::Vector3d x = igl::quadprog<double,3>(H,f,lb,ub);
    REQUIRE(x(0)==0.0);
    REQUIRE(x(1)==0.5);
    REQUIRE(x(1)==0.5);
  }
  {
    Eigen::Matrix3d H = (Eigen::Matrix3d(3,3)<<
        1.06020279605748, 0.387347953430924,-0.653587847224834,
        0.387347953430924, 0.323001970642631, -0.80259721932688,
        -0.653587847224834, -0.80259721932688,  2.73709523329989
        ).finished();
    Eigen::Vector3d f = (Eigen::Vector3d(3,1)<<
        -0.0155804986732503,
        -0.00161174173383921,
        -0.00903647945917485
        ).finished();
    Eigen::Vector3d lb = (Eigen::Vector3d(3,1)<<
        0,
        0,
        0
        ).finished();
    Eigen::Vector3d ub = (Eigen::Vector3d(3,1)<<
        0.015625,
        0.015625,
        0.015625
        ).finished();
	// Windows needs template args spelled out
    Eigen::Vector3d x = igl::quadprog<double,3>(H,f,lb,ub);
    Eigen::Vector3d xexact(0.015625, 0.013732474124087, 0.0110593284260843);
    REQUIRE((x-xexact).array().abs().maxCoeff() < 1e-4);
  }
}

TEST_CASE("quadprog: box2", "[igl]" )
{
  {
    Eigen::Matrix2d H = (Eigen::Matrix2d(2,2)<<
        0.683698654982294,-0.0521997092763332,
        -0.0521997092763332,  0.800535063458999
        ).finished();
    Eigen::Vector2d f = (Eigen::Vector2d(2,1)<<
        -0.733716332234936,
        2.56401312736278
        ).finished();
    Eigen::Vector2d lb = (Eigen::Vector2d(2,1)<<
        -1,
        -1
        ).finished();
    Eigen::Vector2d ub = (Eigen::Vector2d(2,1)<<
        1,
        1
        ).finished();
	// Windows needs template args spelled out
    Eigen::Vector2d x = igl::quadprog<double,2>(H,f,lb,ub);
    //std::cout<<igl::matlab_format(x,"x")<<std::endl;
    REQUIRE(abs(x(0)-0.99680848864073357)<1e-7);
    REQUIRE(abs(x(1)- -1.)<1e-7);
  }
  {
    Eigen::Matrix2d H = (Eigen::Matrix2d(2,2)<<
        0.0708025678527926,-0.158030756288795,
        -0.158030756288795, 0.360468620664163
        ).finished();
    Eigen::Vector2d f = (Eigen::Vector2d(2,1)<<
        0.207229875768196,
        -0.547595351878845
        ).finished();
    Eigen::Vector2d lb = (Eigen::Vector2d(2,1)<<
        -1,
        -1
        ).finished();
    Eigen::Vector2d ub = (Eigen::Vector2d(2,1)<<
        1,
        1
        ).finished();
	// Windows needs template args spelled out
    Eigen::Vector2d x = igl::quadprog<double,2>(H,f,lb,ub);
    //std::cout<<igl::matlab_format(x,"x")<<std::endl;
    REQUIRE(abs(x(0)- -0.69487761491492939)<1e-7);
    REQUIRE(abs(x(1)-1.0)<1e-7);
  }
  {
    Eigen::Matrix2d H;
    H<<0.4000,-0.2000,-0.2000,1.0000;
    Eigen::Vector2d f(-0.3000,4.0000);
    Eigen::Vector2d lb(0,0);
    Eigen::Vector2d ub(1,1);
	// Windows needs template args spelled out
    Eigen::Vector2d x = igl::quadprog<double,2>(H,f,lb,ub);
    //std::cout<<igl::matlab_format(x,"x")<<std::endl;
    REQUIRE(abs(x(0)-0.75)<2e-16);
    REQUIRE(abs(x(1)-0.0)<2e-16);
  }
}

