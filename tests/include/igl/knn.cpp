#include <test_common.h>
#include <igl/list_to_matrix.h>
#include <igl/knn.h>
#include <igl/octree.h>



TEST_CASE("knn", "[igl]")
{

    auto run_knn = [](int k, const Eigen::MatrixXd & P, const Eigen::MatrixXd& V, typename Eigen::MatrixXi& I) {
        std::vector<std::vector<int> > point_indices;
        Eigen::Matrix<int,Eigen::Dynamic,8> CH;
        Eigen::Matrix<double,Eigen::Dynamic,3> CN;
        Eigen::Matrix<double,Eigen::Dynamic,1> W;
        Eigen::Matrix<double,Eigen::Dynamic,1> A;

        igl::octree(V,point_indices,CH,CN,W);

        igl::knn(P,V,k,point_indices,CH,CN,W,I);
    };
    Eigen::Matrix<double,Eigen::Dynamic,3> V;
    Eigen::Matrix<double,Eigen::Dynamic,3> P;
    Eigen::MatrixXi I;
    Eigen::MatrixXi answer;
    {
        //Test some simple points off a unit cube
        V.resize(8,3);
        V << 
            0,0,1,
            0,1,0,
            0,1,1,
            0,0,0,
            1,0,0,
            1,1,0,
            1,1,1,
            1,0,1;

        P.resize(3,3);

        P << 
            0 ,0 ,0.6,
              0.3 ,0.1 ,0.2,
              .7,.6,0;


        answer.resize(3,3);
        answer << 0,3,2,
               3,4,0,
               5,4,1;
        run_knn(3,P,V,I);
        REQUIRE (I == answer);

    }
    {
        //Test whether the kdtree works when passed things of different size

        V.resize(2,3);
        V << 0,0,0,
          1,1,1;

        P << 0,0,0,
            -1,-1,-1,
            2,2,2;

        run_knn(10,P,V,I);
        answer.resize(3,2);
        answer << 
            0,1,
            0,1,
            1,0;

        REQUIRE (I == answer);
    }



}
