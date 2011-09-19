//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 20 Sept 2011

#ifndef IGL_WRITEOBJ_H
#define IGL_WRITEOBJ_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   V  eigen double matrix #V by 3 (mesh vertices)
  //   F  eigen int matrix #F by 3 (mesh indices)
  // Returns true on success, false on error
  bool writeOBJ(const std::string str, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
}

// Implementation
#include <iostream>
#include <fstream>
#include <cstdio>

bool igl::writeOBJ(const std::string str, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
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
#endif
