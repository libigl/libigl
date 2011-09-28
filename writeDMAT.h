#ifndef IGL_WRITEDMAT_H
#define IGL_WRITEDMAT_H
// See writeDMAT.h for a description of the .dmat file type
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Write a matrix using ascii dmat file type
  //
  // Inputs:
  //   file_name  path to .dmat file
  //   W  eigen matrix containing to-be-written coefficients
  // Returns true on success, false on error
  //
  inline bool writeDMAT(const std::string file_name, const Eigen::MatrixXd & W);
}

// Implementation
#include <cstdio>

inline bool igl::writeDMAT(const std::string file_name, const Eigen::MatrixXd & W)
{
  FILE * fp = fopen(file_name.c_str(),"w");
  if(fp == NULL)
  {
    fclose(fp);
    fprintf(stderr,"IOError: writeDMAT() could not open %s...",file_name.c_str());
    return false; 
  }
  // first line contains number of rows and number of columns
  fprintf(fp,"%d %d\n",(int)W.cols(),(int)W.rows());
  // Loop over columns slowly
  for(int j = 0;j < W.cols();j++)
  {
    // loop over rows (down columns) quickly
    for(int i = 0;i < W.rows();i++)
    {
      fprintf(fp,"%lg\n",W(i,j));
    }
  }
  fclose(fp);
  return true;
}
#endif
