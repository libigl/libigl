#ifndef IGL_WRITEDMAT_H
#define IGL_WRITEDMAT_H
// See writeDMAT.h for a description of the .dmat file type
#include <string>
namespace igl
{
  // Write a matrix using ascii dmat file type
  //
  // Template:
  //   Mat  matrix type that supports .rows(), .cols(), operator(i,j)
  // Inputs:
  //   file_name  path to .dmat file
  //   W  eigen matrix containing to-be-written coefficients
  // Returns true on success, false on error
  //
  template <class Mat>
  inline bool writeDMAT(const std::string file_name, const Mat & W);
}

// Implementation
#include <cstdio>

  template <class Mat>
inline bool igl::writeDMAT(const std::string file_name, const Mat & W)
{
  FILE * fp = fopen(file_name.c_str(),"w");
  if(fp == NULL)
  {
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
      fprintf(fp,"%lg\n",(double)W(i,j));
    }
  }
  fclose(fp);
  return true;
}
#endif
