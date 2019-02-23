#include <test_common.h>
#include <igl/unique_simplices.h>

TEST_CASE("igl_unique_simplices: duplicate_triangles", "[igl]")
{
  const Eigen::MatrixXi F = (Eigen::MatrixXi(2,3)<<0,1,2,0,1,2).finished();

  // All possible permutations of the same triangle
  for(int di = -1;di<2;di+=2)
  {
    for(int dj = -1;dj<2;dj+=2)
    {
      for(int i = 0;i<3;i++)
      {
        for(int j = 0;j<3;j++)
        {
          Eigen::MatrixXi Fij = F;
          for(int c = 0;c<3;c++)
          {
            Fij(0,c) = (Fij(0,c)+3+di*i)%3;
            Fij(1,c) = (Fij(1,c)+3+dj*j)%3;
          }
          Eigen::MatrixXi Fu;
          Eigen::VectorXi IA,IC;
          igl::unique_simplices(Fij,Fu,IA,IC);
          // There's only one unique simplex
          REQUIRE (1 == Fu.rows());
        }
      }
    }
  }

}
