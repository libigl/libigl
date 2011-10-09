#ifndef IGL_READDMAT_H
#define IGL_READDMAT_H
// .dmat is a simple ascii matrix file type, defined as follows. The first line
// is always:
// <#columns> <#rows>
// Then the coefficients of the matrix are given separated by whitespace with
// columns running fastest.
//
// Example:
//   The matrix m = [1 2 3; 4 5 6];
//   corresponds to a .dmat file containing:
//   3 2
//   1 4 2 5 3 6
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Read a matrix from an ascii dmat file
  //
  // Inputs:
  //   file_name  path to .dmat file
  // Outputs:
  //   W  eigen matrix containing read-in coefficients
  // Returns true on success, false on error
  //
  inline bool readDMAT(const std::string file_name, Eigen::MatrixXd & W);
}

// Implementation
#include "verbose.h"
#include <cstdio>

inline bool igl::readDMAT(const std::string file_name, Eigen::MatrixXd & W)
{
  FILE * fp = fopen(file_name.c_str(),"r");
  if(fp == NULL)
  {
    fprintf(stderr,"IOError: readDMAT() could not open %s...\n",file_name.c_str());
    return false; 
  }
  // first line contains number of rows and number of columns
  int num_cols, num_rows;
  int res = fscanf(fp,"%d %d\n",&num_cols,&num_rows);
  if(res != 2)
  {
    fclose(fp);
    fprintf(stderr,"IOError: readDMAT() first row should be [num cols] [num rows]...\n");
    return false;
  }

  // check that number of columns and rows are sane
  if(num_cols < 0)
  {
    fclose(fp);
    fprintf(stderr,"IOError: readDMAT() number of columns %d < 0\n",num_cols);
    return false;
  }
  if(num_rows < 0)
  {
    fclose(fp);
    fprintf(stderr,"IOError: readDMAT() number of rows %d < 0\n",num_rows);
    return false;
  }

  // Resize output to fit matrix
  W.resize(num_rows,num_cols);

  // Loop over columns slowly
  for(int j = 0;j < num_cols;j++)
  {
    // loop over rows (down columns) quickly
    for(int i = 0;i < num_rows;i++)
    {
      if(fscanf(fp," %lg",&W(i,j)) != 1)
      {
        fclose(fp);
        fprintf(
          stderr,
          "IOError: readDMAT() bad format after reading %d entries\n",
          j*num_rows + i);
        return false;
      }
    }
  }
  fclose(fp);
  return true;
}
#endif
