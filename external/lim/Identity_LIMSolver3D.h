// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef IDENTITY_LIM_SOLVER_3D_H
#define IDENTITY_LIM_SOLVER_3D_H

#include "LIMSolver3D.h"

class Identity_LIMSolver3D : public LIMSolver3D
{
public:

  Identity_LIMSolver3D();
  ~Identity_LIMSolver3D();

private:
  
  Eigen::SparseMatrix<double> Identity;

  void debugOutput(std::stringstream& info);
  void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx);
  double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess);
};

#endif