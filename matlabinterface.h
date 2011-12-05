//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// WARNING: These functions require matlab installed
// Additional header folder required:
//   /Applications/MATLAB_R2011a.app/extern/include
// Additional binary lib to be linked with:
// /Applications/MATLAB_R2011a.app/bin/maci64/libeng.dylib
// /Applications/MATLAB_R2011a.app/bin/maci64/libmx.dylib

// MAC ONLY:
// Add to the environment variables:
// DYLD_LIBRARY_PATH = /Applications/MATLAB_R2011a.app/bin/maci64
// PATH = /opt/local/bin:/opt/local/sbin:/Applications/MATLAB_R2011a.app/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/texbin:/usr/X11/bin

#ifndef MATLAB_INTERFACE_H
#define MATLAB_INTERFACE_H

#include <Eigen/Core>
#include <string>

#include <complex>
#include <cassert>
#include <map>
#include <string>
#include <vector>

#include "engine.h"  // Matlab engine header

namespace igl 
{
  // Init the MATLAB engine 
  // (no need to call it directly since it is automatically invoked by any other command)
  inline void mlinit(Engine** engine);
  
  // Closes the MATLAB engine
  inline void mlclose(Engine** engine);
  
  // Send a matrix to MATLAB
  inline void mlsetmatrix(Engine** engine, std::string name, const Eigen::MatrixXd& M);
  
  // Send a matrix to MATLAB
  inline void mlsetmatrix(Engine** engine, std::string name, const Eigen::MatrixXf& M);

  // Send a matrix to MATLAB
  inline void mlsetmatrix(Engine** engine, std::string name, const Eigen::MatrixXi& M);
  
  // Receive a matrix from MATLAB
  inline void mlgetmatrix(Engine** engine, std::string name, Eigen::MatrixXd& M);
  
  // Receive a matrix from MATLAB
  inline void mlgetmatrix(Engine** engine, std::string name, Eigen::MatrixXf& M);

  // Receive a matrix from MATLAB
  inline void mlgetmatrix(Engine** engine, std::string name, Eigen::MatrixXi& M);
  
  // Send a single scalar to MATLAB
  inline void mlsetscalar(Engine** engine, std::string name, double s);
  
  // Receive a single scalar from MATLAB
  inline double mlgetscalar(Engine** engine, std::string name);
  
  // Execute arbitrary MATLAB code and return the MATLAB output
  inline std::string mleval(Engine** engine, std::string code);
  
}

//Implementation

// Init the MATLAB engine 
// (no need to call it directly since it is automatically invoked by any other command)
inline void igl::mlinit(Engine** mlengine)
{
  *mlengine = engOpen("\0");
}

// Closes the MATLAB engine
inline void igl::mlclose(Engine** mlengine)
{
  engClose(*mlengine);
  *mlengine = 0;
}

// Send a matrix to MATLAB
inline void igl::mlsetmatrix(Engine** mlengine, std::string name, const Eigen::MatrixXd& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  mxArray *A = mxCreateDoubleMatrix(M.rows(), M.cols(), mxREAL);
  double *pM = mxGetPr(A);
  
  int c = 0;
  for(int j=0; j<M.cols();++j)
    for(int i=0; i<M.rows();++i)
      pM[c++] = double(M(i,j));
  
  engPutVariable(*mlengine, name.c_str(), A);
  mxDestroyArray(A);
}

// Send a matrix to MATLAB
inline void igl::mlsetmatrix(Engine** mlengine, std::string name, const Eigen::MatrixXf& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  mxArray *A = mxCreateDoubleMatrix(M.rows(), M.cols(), mxREAL);
  double *pM = mxGetPr(A);
  
  int c = 0;
  for(int j=0; j<M.cols();++j)
    for(int i=0; i<M.rows();++i)
      pM[c++] = double(M(i,j));
  
  engPutVariable(*mlengine, name.c_str(), A);
  mxDestroyArray(A);
}

// Send a matrix to MATLAB
inline void igl::mlsetmatrix(Engine** mlengine, std::string name, const Eigen::MatrixXi& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  mxArray *A = mxCreateDoubleMatrix(M.rows(), M.cols(), mxREAL);
  double *pM = mxGetPr(A);
  
  int c = 0;
  for(int j=0; j<M.cols();++j)
    for(int i=0; i<M.rows();++i)
      pM[c++] = double(M(i,j))+1;
  
  engPutVariable(*mlengine, name.c_str(), A);
  mxDestroyArray(A);
}

// Receive a matrix from MATLAB
inline void igl::mlgetmatrix(Engine** mlengine, std::string name, Eigen::MatrixXd& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  unsigned long m = 0;
  unsigned long n = 0;
  std::vector<double> t;
  
  mxArray *ary = engGetVariable(*mlengine, name.c_str());
  if (ary == NULL)
  {
    m = 0;
    n = 0;
    M = Eigen::MatrixXd(0,0);
  }
  else
  {
    m = mxGetM(ary);
    n = mxGetN(ary);
    M = Eigen::MatrixXd(m,n);
    
    double *pM = mxGetPr(ary);
    
    int c = 0;
    for(int j=0; j<M.cols();++j)
      for(int i=0; i<M.rows();++i)
        M(i,j) = pM[c++];
  }
  
  mxDestroyArray(ary);
}

inline void igl::mlgetmatrix(Engine** mlengine, std::string name, Eigen::MatrixXf& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  unsigned long m = 0;
  unsigned long n = 0;
  std::vector<double> t;
  
  mxArray *ary = engGetVariable(*mlengine, name.c_str());
  if (ary == NULL)
  {
    m = 0;
    n = 0;
    M = Eigen::MatrixXf(0,0);
  }
  else
  {
    m = mxGetM(ary);
    n = mxGetN(ary);
    M = Eigen::MatrixXf(m,n);
    
    double *pM = mxGetPr(ary);
    
    int c = 0;
    for(int j=0; j<M.cols();++j)
      for(int i=0; i<M.rows();++i)
        M(i,j) = pM[c++];
  }
  
  mxDestroyArray(ary);
}

// Receive a matrix from MATLAB
inline void igl::mlgetmatrix(Engine** mlengine, std::string name, Eigen::MatrixXi& M)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  unsigned long m = 0;
  unsigned long n = 0;
  std::vector<double> t;
  
  mxArray *ary = engGetVariable(*mlengine, name.c_str());
  if (ary == NULL)
  {
    m = 0;
    n = 0;
    M = Eigen::MatrixXi(0,0);
  }
  else
  {
    m = mxGetM(ary);
    n = mxGetN(ary);
    M = Eigen::MatrixXi(m,n);
    
    double *pM = mxGetPr(ary);
    
    int c = 0;
    for(int j=0; j<M.cols();++j)
      for(int i=0; i<M.rows();++i)
        M(i,j) = int(pM[c++])-1;
  }
  
  mxDestroyArray(ary);
}

// Send a single scalar to MATLAB
inline void igl::mlsetscalar(Engine** mlengine, std::string name, double s)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  Eigen::MatrixXd M(1,1);
  M(0,0) = s;
  mlsetmatrix(mlengine, name, M);
}

// Receive a single scalar from MATLAB
inline double igl::mlgetscalar(Engine** mlengine, std::string name)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  Eigen::MatrixXd M;
  mlgetmatrix(mlengine, name,M);
  return M(0,0);
}

// Execute arbitrary MATLAB code and return the MATLAB output
inline std::string igl::mleval(Engine** mlengine, std::string code)
{
  if (*mlengine == 0)
    mlinit(mlengine);
  
  const char *matlab_code = code.c_str();
  const int BUF_SIZE = 4096*4096;
  // allocate on the heap to avoid running out of stack
  std::string bufauto(BUF_SIZE+1, '\0');
  char *buf = &bufauto[0];
  
  assert(matlab_code != NULL);
  
  // Use RAII ensure that on leaving this scope, the output buffer is
  // always nullified (to prevent Matlab from accessing memory that might
  // have already been deallocated).
  struct cleanup {
    Engine *m_ep;
    cleanup(Engine *ep) : m_ep(ep) { }
    ~cleanup() { engOutputBuffer(m_ep, NULL, 0); }
  } cleanup_obj(*mlengine);
  
  if (buf != NULL)
    engOutputBuffer(*mlengine, buf, BUF_SIZE);
  
  int res = engEvalString(*mlengine, matlab_code);
  
  if (res != 0) {
    std::ostringstream oss;
    oss << "ERROR: Matlab command failed with error code " << res << ".\n";
    return oss.str();
  }
  
  if (buf[0] == '>' && buf[1] == '>' && buf[2] == ' ')
    buf += 3;
  if (buf[0] == '\n') ++buf;
  
  return std::string(buf);
}


#endif