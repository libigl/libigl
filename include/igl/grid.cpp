#include "grid.h"

IGL_INLINE void igl::grid(const Eigen::RowVector3i & res, Eigen::MatrixXd & GV)
{
  using namespace Eigen;
  GV.resize(res(0)*res(1)*res(2),3);
  for(int zi = 0;zi<res(2);zi++)
  {
    const auto lerp = 
      [&](const double di, const int d)->double{return di/(double)(res(d)-1);};
    const double z = lerp(zi,2);
    for(int yi = 0;yi<res(1);yi++)
    {
      const double y = lerp(yi,1);
      for(int xi = 0;xi<res(0);xi++)
      {
        const double x = lerp(xi,0);
        GV.row(xi+res(0)*(yi + res(1)*zi)) = RowVector3d(x,y,z);
      }
    }
  }
}

