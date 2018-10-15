#include <test_common.h>

#include <CGAL/Simple_cartesian.h>
#include <igl/copyleft/cgal/hausdorff.h>
#include <igl/copyleft/cgal/point_mesh_squared_distance.h>
#include <igl/upsample.h>

TEST(hausdorff, knightVScheburashka) 
{
  Eigen::MatrixXd VA,VB;
  Eigen::MatrixXi FA,FB;
  test_common::load_mesh("decimated-knight.obj", VA, FA);
  test_common::load_mesh("cheburashka.off", VB, FB);
  //typedef CGAL::Epeck Kernel;
  typedef CGAL::Simple_cartesian<double> Kernel;
  CGAL::AABB_tree<
    CGAL::AABB_traits<Kernel, 
      CGAL::AABB_triangle_primitive<Kernel, 
        typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
      >
    >
  > treeB;
  std::vector<CGAL::Triangle_3<Kernel> > TB;
  {
    igl::copyleft::cgal::point_mesh_squared_distance_precompute(VB,FB,treeB,TB);
  }
  std::vector<Eigen::VectorXd> U(5);
  for(int j = 0;j<U.size();j++)
  {
    if(j>0)
    {
      igl::upsample(Eigen::MatrixXd(VA),Eigen::MatrixXi(FA),VA,FA);
    }
    const int m = FA.rows();
    U[j].resize(m);
    for(int i = 0;i<m;i++)
    {
      Eigen::MatrixXd Vi(3,3);
      Vi<<VA.row(FA(i,0)),VA.row(FA(i,1)),VA.row(FA(i,2));
      double l;
      igl::copyleft::cgal::hausdorff(Vi,treeB,TB,l,U[j](i));
      if(j>0&&i%4==3)
      {
        double u4 = std::numeric_limits<double>::infinity();
        for(int u = 0;u<4;u++)
        {
          u4 = std::min(u4,U[j](i-u));
        }
        ASSERT_LE(u4,U[j-1](i/4));
      }
    }
    break;
  }
}

