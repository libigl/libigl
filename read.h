//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011


#ifndef IGL_READ_H
#define IGL_READ_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
    // read mesh from an ascii file with automatic detection of file format. supported: obj, off)
  // Inputs:
  //   str  path to .obj/.off file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  bool read(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}

// Implementation
#include <readOBJ.h>
#include <readOFF.h>
bool igl::read(const std::string str, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    const char* p;
    for (p = str.c_str(); *p != '\0'; p++)
        ;
    while (*p != '.')
        p--;
    
    if (!strcmp(p, ".obj") || !strcmp(p, ".OBJ"))
    {
        return igl::readOBJ(str,V,F);
    }else if (!strcmp(p, ".off") || !strcmp(p, ".OFF"))
    {
        return igl::readOFF(str,V,F);
    }else
    {
      fprintf(stderr,"read() does not recognize extension: %s\n",p);
      return false;
    }
}

#endif
