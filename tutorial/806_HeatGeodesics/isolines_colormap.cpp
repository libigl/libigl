#include "isolines_colormap.h"
#include <igl/isolines_map.h>


Eigen::MatrixXd isolines_colormap()
{
  const int num_intervals = 30;
  Eigen::MatrixXd CM(num_intervals,3);
  // Colormap texture
  for(int i = 0;i<num_intervals;i++)
  {
    double t = double(num_intervals - i - 1)/double(num_intervals-1);
    CM(i,0) = std::max(std::min(2.0*t-0.0,1.0),0.0);
    CM(i,1) = std::max(std::min(2.0*t-1.0,1.0),0.0);
    CM(i,2) = std::max(std::min(6.0*t-5.0,1.0),0.0);
  }
  igl::isolines_map(Eigen::MatrixXd(CM),CM);
  return CM;
}
