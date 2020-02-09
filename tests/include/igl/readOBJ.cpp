#include <igl/readOBJ.h>
#include <test_common.h>
#include <iostream>
#include <string>
#include <tuple>

TEST_CASE("readOBJ: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOBJ
    // directly...
    igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (F.rows() == 12);
}

TEST_CASE("readOBJ: Obj with material", "[igl]")
{
    std::vector<std::vector<double > > V;
    std::vector<std::vector<double > > TC;
    std::vector<std::vector<double > > N;
    std::vector<std::vector<int > > F;
    std::vector<std::vector<int > > FTC;
    std::vector<std::vector<int > >  FN;
    std::vector<std::tuple<std::string, int, int>> FM;
    igl::readOBJ(test_common::data_path("cubewithmaterial.obj"), V, TC, N, F, FTC, FN, FM);

    REQUIRE (V.size() == 8);
    REQUIRE (F.size() == 6);
    for ( const auto& i : FM ) {
        std::cout << "material ";
        std::cout << std::get<0>(i) << ' ';
        std::cout << "fstart ";
        std::cout << std::get<1>(i) << ' ';
        std::cout << "fend ";
        std::cout << std::get<2>(i) << ' ';
        std::cout << std::endl;
    }
    REQUIRE (FM.size() == 2);
}
