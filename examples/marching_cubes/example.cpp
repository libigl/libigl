#include <igl/marching_cubes.h>
#include <iostream>
#include <cstdio>

#include "igl/MCTables.hh"

typedef float ScalarType;
typedef unsigned IndexType;
int main(int argc, char * argv[])
{
  Eigen::Matrix<ScalarType, 1, 3> bb_min, bb_max;
  bb_min<<0.,0.,0.;
  bb_max<<1.,1.,1.;
  
  //diagonal
  Eigen::Matrix<ScalarType, 1, 3> diff = bb_max - bb_min;
  
  
  //center of the sphere
  Eigen::Matrix<ScalarType, 1, 3> center = 0.5 * (bb_min + bb_max);
  
  IndexType xres, yres, zres;
  xres = yres = zres = 10;
  ScalarType radius = 0.42;

  //steps in x,y,z direction
  ScalarType dx = diff[0] / (ScalarType)(xres-1);
  ScalarType dy = diff[1] / (ScalarType)(yres-1);
  ScalarType dz = diff[2] / (ScalarType)(zres-1);
  
  
  Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> points(xres*yres*zres,3);
  Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> values(xres*yres*zres,1);
  
  Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> vertices;
  
  Eigen::Matrix<IndexType, Eigen::Dynamic, 3> faces;

  std::cerr<<"Sphere -- construct grid"<<std::endl;
  for (unsigned int x=0; x<xres; ++x)
    for (unsigned int y=0; y<yres; ++y)
      for (unsigned int z=0; z<zres; ++z)
      {
        int index = x + y*xres + z*xres*yres;
        points.row(index) = bb_min + 
        ScalarType(x)*Eigen::Matrix<ScalarType, 1, 3>(dx,0.,0.) + 
        ScalarType(y)*Eigen::Matrix<ScalarType, 1, 3>(0.,dy,0.) + 
        ScalarType(z)*Eigen::Matrix<ScalarType, 1, 3>(0.,0.,dz);
        
        values[index] = (points.row(index) - center).squaredNorm() - radius*radius;
      }
  
  
  std::cerr<<"Sphere -- marching cubes"<<std::endl;
  igl::marching_cubes(values, 
                      points, 
                      xres, 
                      yres, 
                      zres, 
                      vertices,
                      faces);

  std::cerr<<"Sphere -- saving"<<std::endl;
  FILE * fid = fopen("sphere.obj","w");
  for (unsigned i = 0; i<vertices.rows(); ++i)
    fprintf(fid,"v %.10g %.10g %.10g\n",vertices(i,0),vertices(i,1),vertices(i,2));
  for (unsigned i = 0; i<faces.rows(); ++i)
    fprintf(fid,"f %d %d %d\n",faces(i,0)+1,faces(i,1)+1,faces(i,2)+1);
  fclose(fid);

  std::cerr<<"Sphere -- done."<<std::endl;
  return 0;
}
