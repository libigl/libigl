// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef LS_Conformal_LIM_SOLVER_2D_H
#define LS_Conformal_LIM_SOLVER_2D_H

#include "LIMSolver2D.h"

class LSConformal_LIMSolver2D : public LIMSolver2D
{
public:

  LSConformal_LIMSolver2D();
  ~LSConformal_LIMSolver2D();

private:
  
  Eigen::SparseMatrix<double> L, A;

  // implementation of abstract functions of LIMSolver2D
  void debugOutput(std::stringstream& info);
  void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx);
  double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess);
};

#endif