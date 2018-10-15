#include <test_common.h>
#include <igl/seam_edges.h>
#include <igl/readOBJ.h>

TEST(seam_edges, tet)
{
  Eigen::MatrixXd V,TC,CN;
  Eigen::MatrixXi F,FTC,FN;
  // Load example mesh: GetParam() will be name of mesh file
  igl::readOBJ(test_common::data_path("tet.obj"), V, TC,CN,F,FTC,FN);
  Eigen::MatrixXi seams,boundaries,foldovers;
  igl::seam_edges(V,TC,F,FTC,seams,boundaries,foldovers);

  Eigen::MatrixXi seams_gt(3,4);
  seams_gt<<
    0,0,1,2,
    3,0,0,2,
    1,0,3,2;
  test_common::assert_eq(seams,seams_gt);
  ASSERT_EQ(boundaries.size(),0);
  ASSERT_EQ(foldovers.size(),0);
}

