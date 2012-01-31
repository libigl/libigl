#include "readDMAT.h"

#include "verbose.h"
#include <cstdio>

IGL_INLINE bool igl::readDMAT(const std::string file_name, Eigen::MatrixXd & W)
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
