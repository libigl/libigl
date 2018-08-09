#include <igl/doublearea.h>
#include <igl/internal_angles.h>
#include <igl/is_irregular_vertex.h>
#include <igl/readOBJ.h>
#include <igl/PI.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  MatrixXd V;
  MatrixXi F;

  igl::readOBJ(TUTORIAL_SHARED_PATH "/horse_quad.obj",V,F);

  // Count the number of irregular vertices, the border is ignored
  vector<bool> irregular = igl::is_irregular_vertex(V,F);

  int vertex_count = V.rows();
  int irregular_vertex_count = 
    std::count(irregular.begin(),irregular.end(),true);
  double irregular_ratio = double(irregular_vertex_count)/vertex_count;

  printf("Irregular vertices: \n%d/%d (%.2f%%)\n",
    irregular_vertex_count,vertex_count, irregular_ratio*100);

  // Compute areas, min, max and standard deviation
  VectorXd area;
  igl::doublearea(V,F,area);
  area = area.array() / 2;

  double area_avg   = area.mean();
  double area_min   = area.minCoeff() / area_avg;
  double area_max   = area.maxCoeff() / area_avg;
  double area_sigma = sqrt(((area.array()-area_avg)/area_avg).square().mean());

  printf("Areas (Min/Max)/Avg_Area Sigma: \n%.2f/%.2f (%.2f)\n",
    area_min,area_max,area_sigma);

  // Compute per face angles, min, max and standard deviation
  MatrixXd angles;
  igl::internal_angles(V,F,angles);
  angles = 360.0 * (angles/(2*igl::PI)); // Convert to degrees

  double angle_avg   = angles.mean();
  double angle_min   = angles.minCoeff();
  double angle_max   = angles.maxCoeff();
  double angle_sigma = sqrt( (angles.array()-angle_avg).square().mean() );

  printf("Angles in degrees (Min/Max) Sigma: \n%.2f/%.2f (%.2f)\n",
    angle_min,angle_max,angle_sigma);

}
