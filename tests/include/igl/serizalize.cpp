#include <test_common.h>
#include <igl/serialize.h>


TEST_CASE("serialize: serialize data", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("truck.obj", V,F);

    REQUIRE(igl::serialize(V, "V", test_common::data_path("truck.serialized")));
    REQUIRE(igl::serialize(F, "F", test_common::data_path("truck.serialized")));

    REQUIRE(V.rows() == 2956);
    REQUIRE(F.rows() == 4770);
}

TEST_CASE("serialize: deserialize data", "[igl]")
{
    Eigen::MatrixXd V, V_orig;
    Eigen::MatrixXi F, F_orig;
    test_common::load_mesh("truck.obj", V_orig,F_orig);
    REQUIRE(V_orig.rows() == 2956);
    REQUIRE(F_orig.rows() == 4770);

    REQUIRE(igl::deserialize(V, "V", test_common::data_path("truck.serialized")));
    REQUIRE(igl::deserialize(F, "F", test_common::data_path("truck.serialized")));

    test_common::assert_eq(V, V_orig);
    test_common::assert_eq(F, F_orig);
}

TEST_CASE("serialize: serialize string", "[igl]")
{
    std::string description = "libigl - A simple C++ geometry processing library.";

    REQUIRE(igl::serialize(description, "description", test_common::data_path("string.serialized")));
}

TEST_CASE("serialize: deserialize string", "[igl]")
{
    std::string deserialized;
    std::string description = "libigl - A simple C++ geometry processing library.";

    REQUIRE(igl::deserialize(deserialized, "description", test_common::data_path("string.serialized")));

    REQUIRE(deserialized == description);
}

TEST_CASE("serialize: cross platform deserialize", "[igl]")
{

  std::string str_stored = "libigl - A simple C++ geometry processing library.";
  std::vector<std::string> str_vec_stored = {"libigl"," - A simple", "C++", "geometry processing library."};

  std::string str_load;
  std::vector<std::string> str_vec_load;

  // load serialized file under ubuntu (gcc 7.4.0)
  igl::TypeCheck option = igl::TypeCheck::ERROR_ON_MISMATCH;
  REQUIRE(igl::deserialize(str_load    ,"str"    ,test_common::data_path("cross_platform_mac.serialized"),option));
  REQUIRE(igl::deserialize(str_vec_load,"str_vec",test_common::data_path("cross_platform_mac.serialized"),option));

  REQUIRE(igl::deserialize(str_load    ,"str"    ,test_common::data_path("cross_platform_windows.serialized"),option));
  REQUIRE(igl::deserialize(str_vec_load,"str_vec",test_common::data_path("cross_platform_windows.serialized"),option));

  //REQUIRE(igl::deserialize(str_ubuntu,"str",test_common::data_path("cross_platform_ubuntu.serialized"),option));
  //REQUIRE(igl::deserialize(str_vec_ubuntu,"str_vec",test_common::data_path("cross_platform_ubuntu.serialized"),option));
  // REQUIRE(igl::deserialize(str_hpc,"str",test_common::data_path("cross_platform_hpc.serialized"),option));
  // REQUIRE(igl::deserialize(str_vec_hpc,"str_vec",test_common::data_path("cross_platform_hpc.serialized"),option));

  // REQUIRE(str_vec_mac == str_vec_stored);
  // REQUIRE(str_vec_windows == str_vec_stored);
}

