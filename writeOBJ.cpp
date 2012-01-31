#include "writeOBJ.h"

#include <iostream>
#include <fstream>
#include <cstdio>

IGL_INLINE bool igl::writeOBJ(const std::string str, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
  std::ofstream s(str.c_str());

  if(!s.is_open())
  {
    fprintf(stderr,"IOError: writeOBJ() could not open %s\n",str.c_str());
    return false;
  }

  for(int i=0;i<V.rows();++i)
      s << "v " << V(i,0) << " " << V(i,1) << " " << V(i,2) << std::endl;
  
  for(int i=0;i<F.rows();++i)
      s << "f " << F(i,0)+1 << " " << F(i,1)+1 << " " << F(i,2)+1 << std::endl;
  
  s.close();
  return true;
}
