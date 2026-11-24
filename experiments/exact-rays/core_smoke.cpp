#include <igl/eytzinger_aabb_ray_intersection.h>
int main()
{
  const Eigen::RowVector3d origin(0,0,-1),direction(0,0,1);
  Eigen::Matrix<double,Eigen::Dynamic,3> low(1,3),high(1,3);
  low<<-1,-1,0; high<<1,1,0;
  Eigen::VectorXi leaf(1); leaf<<0;
  std::vector<int> hits;
  igl::eytzinger_aabb_ray_intersection<true>(origin,direction,
    [](int){return true;},[](int,int){return false;},low,high,leaf,hits);
  return hits!=std::vector<int>{0};
}
