#include <test_common.h>

TEST_CASE("readOBJ: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOBJ
    // directly...
    test_common::load_mesh("cube.obj", V, F);
    REQUIRE (V.rows() == 8);
    REQUIRE (F.rows() == 12);
}

TEST_CASE("readOBJ: Obj with material", "[igl]")
{
    //Eigen::MatrixXd V;
    //Eigen::MatrixXi F;
    // wait... so this is actually testing test_common::load_mesh not readOBJ
    // directly...

    std::vector<std::vector<double > > V;
    std::vector<std::vector<double > > TC;
    std::vector<std::vector<double > > N;
    std::vector<std::vector<int > > F;
    std::vector<std::vector<int > > FTC;
    std::vector<std::vector<int > >  FN;
    std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::string> > FM;
    test_common::load_obj_with_material<double, int>("cube.obj", V, TC, N, F, FTC, FN, FM);
    REQUIRE (V.size() == 8);
    //REQUIRE (F.rows() == 12);
}
