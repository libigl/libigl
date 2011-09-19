//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#ifndef IGL_READOFF_H
#define IGL_READOFF_H

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // read mesh from a ascii off file
  // Inputs:
  //   str  path to .off file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  bool readOFF (const std::string meshfile, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}

// Implementation

bool igl::readOFF (const std::string meshfile, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    int vnum, fnum;
    FILE *fp = fopen (meshfile.c_str(), "r");
    
    if (!fp)
    {
      fprintf (stderr, "readOFF(): could not open file %s", meshfile.c_str());
      return false;
    }
    
    fscanf (fp, "OFF\n%d %d 0\n",  &vnum, &fnum);
    
    V = Eigen::MatrixXd (vnum, 3);
    F = Eigen::MatrixXi (fnum, 3);
    
    for (unsigned i = 0; i < V.rows(); i++)
        fscanf (fp, "%lf %lf %lf\n", &V(i,0), &V(i,1), &V(i,2));
    
    for (unsigned i = 0; i < F.rows(); i++)
        fscanf (fp, "3 %d %d %d\n", &F(i,0), &F(i,1), &F(i,2));
    
    fclose (fp);
    return true;
}

#endif
