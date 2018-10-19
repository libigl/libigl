#include <test_common.h>
#include <igl/tet_tet_adjacency.h>
#include <igl/readMESH.h>
#include <iostream>


// class tet_tet_adjacency : public ::testing::TestWithParam<std::string> {};

// TEST_P(tet_tet_adjacency, dot)
// {
//   Eigen::MatrixXd V;
//   Eigen::MatrixXi F, T, TT,TTi;
//   // Load example mesh: GetParam() will be name of mesh file
//   igl::readMESH(test_common::data_path(GetParam()), V, T, F);
//   igl::tet_tet_adjacency(T, TT, TTi);
//   REQUIRE (TT.rows() == T.rows());
//   REQUIRE (TTi.rows() == T.rows());
//   REQUIRE (TT.cols() == T.cols());
//   REQUIRE (TTi.cols() == T.cols());
//   for(int t = 0;t<T.rows();t++)
//   {
//     for(int c = 0; c<4 ;c++)
//     {
//       if(TT(t, c) >= 0)
//       {
//         REQUIRE (T.rows() > TT(t, c));
//         REQUIRE (0 <= TTi(t, c));
//         REQUIRE (4 > TTi(t, c));
//         REQUIRE (t == TT(TT(t, c), TTi(t,c)));
//       }
//     }
//   }
// }

// INSTANTIATE_TEST_CASE_P
// (
//   tet_meshes,
//   tet_tet_adjacency,
//   ::testing::ValuesIn(test_common::tet_meshes()),
//   test_common::string_test_name
// );
