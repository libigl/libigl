//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_WRITE_H
#define IGL_WRITE_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // write mesh to an ascii file with automatic detection of file format. supported: obj, off)
  // Known Bugs:
  //  Does not correctly find file extensions: myfile.foo.off 
  inline bool write(
    const std::string str, 
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi& F);
}

// Implementation
#include <writeOBJ.h>
#include <writeOFF.h>

inline bool igl::write(
  const std::string str, 
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& F)
{
  const char* p;
  for (p = str.c_str(); *p != '\0'; p++)
    ;
  while (*p != '.')
    p--;
  
  if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
    return igl::writeOBJ(str,V,F);
  
  if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
    return igl::writeOFF(str,V,F);
}

#endif
