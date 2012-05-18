#include "readDMAT.h"

#include "verbose.h"
#include <cstdio>
#include <iostream>
#include <cassert>

// Static helper method reads the first to elements in the given file
// Inputs:
//   fp  file pointer of .dmat file that was just opened
// Outputs:
//   num_rows  number of rows
//   num_cols number of columns
// Returns 
//   0  success
//   1  did not find header
//   2  bad num_cols
//   3  bad num_rows
static inline int readDMAT_read_header(FILE * fp, int & num_rows, int & num_cols)
{
  // first line contains number of rows and number of columns
  int res = fscanf(fp,"%d %d\n",&num_cols,&num_rows);
  if(res != 2)
  {
    return 1;
  }
  // check that number of columns and rows are sane
  if(num_cols < 0)
  {
    fprintf(stderr,"IOError: readDMAT() number of columns %d < 0\n",num_cols);
    return 2;
  }
  if(num_rows < 0)
  {
    fprintf(stderr,"IOError: readDMAT() number of rows %d < 0\n",num_rows);
    return 3;
  }
  return 0;
}

#ifndef IGL_NO_EIGEN
template <typename DerivedW>
IGL_INLINE bool igl::readDMAT(const std::string file_name,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  FILE * fp = fopen(file_name.c_str(),"r");
  if(fp == NULL)
  {
    fprintf(stderr,"IOError: readDMAT() could not open %s...\n",file_name.c_str());
    return false; 
  }
  int num_rows,num_cols;
  int head_success = readDMAT_read_header(fp,num_rows,num_cols);
  if(head_success != 0)
  {
    if(head_success == 1)
    {
      fprintf(stderr,
        "IOError: readDMAT() first row should be [num cols] [num rows]...\n");
    }
    fclose(fp);
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
      double d;
      if(fscanf(fp," %lg",&d) != 1)
      {
        fclose(fp);
        fprintf(
          stderr,
          "IOError: readDMAT() bad format after reading %d entries\n",
          j*num_rows + i);
        return false;
      }
      W(i,j) = d;
    }
  }

  fclose(fp);
  return true;
}
#endif

template <typename Scalar>
IGL_INLINE bool igl::readDMAT(
  const std::string file_name, 
  std::vector<std::vector<Scalar> > & W)
{
  FILE * fp = fopen(file_name.c_str(),"r");
  if(fp == NULL)
  {
    fprintf(stderr,"IOError: readDMAT() could not open %s...\n",file_name.c_str());
    return false; 
  }
  int num_rows,num_cols;
  bool head_success = readDMAT_read_header(fp,num_rows,num_cols);
  if(head_success != 0)
  {
    if(head_success == 1)
    {
      fprintf(stderr,
        "IOError: readDMAT() first row should be [num cols] [num rows]...\n");
    }
    fclose(fp);
    return false;
  }

  // Resize for output
  W.resize(num_rows,typename std::vector<Scalar>(num_cols));

  // Loop over columns slowly
  for(int j = 0;j < num_cols;j++)
  {
    // loop over rows (down columns) quickly
    for(int i = 0;i < num_rows;i++)
    {
      double d;
      if(fscanf(fp," %lg",&d) != 1)
      {
        fclose(fp);
        fprintf(
          stderr,
          "IOError: readDMAT() bad format after reading %d entries\n",
          j*num_rows + i);
        return false;
      }
      W[i][j] = (Scalar)d;
    }
  }

  // Try to read header for binary part
  head_success = readDMAT_read_header(fp,num_rows,num_cols);
  if(head_success == 0)
  {
    assert(W.size() == 0);
    // Resize for output
    W.resize(num_rows,typename std::vector<Scalar>(num_cols));
    double * Wraw = new double[num_rows*num_cols];
    fread(Wraw, 8, num_cols*num_rows, fp);
    // Loop over columns slowly
    for(int j = 0;j < num_cols;j++)
    {
      // loop over rows (down columns) quickly
      for(int i = 0;i < num_rows;i++)
      {
        W[i][j] = Wraw[j*num_rows+i];
      }
    }
  }

  fclose(fp);
  return true;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template bool igl::readDMAT<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
