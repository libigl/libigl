// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef POISSON_LIM_SOLVER_2D_H
#define POISSON_LIM_SOLVER_2D_H

#include "LIMSolver2D.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

class Poisson_LIMSolver2D : public LIMSolver2D
{
public:

  Poisson_LIMSolver2D();
  ~Poisson_LIMSolver2D();

  // E = 0.5*||Gx-b||^2 = 0.5*x'G'Gx - b'Gx
  Eigen::Matrix<double,Eigen::Dynamic,1> b;

private:

  double constantEnergyPart;
  
  Eigen::SparseMatrix<double> G;
  Eigen::SparseMatrix<double> L;
  Eigen::Matrix<double,Eigen::Dynamic,1> GTb;

  Eigen::Matrix<double*,1,Eigen::Dynamic> lapcianHessianCoeffs;
  
  // implementation of abstract functions of LIMSolver2D
  void debugOutput(std::stringstream& info);
  void getProblemSize();
  void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx);
  double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess);
};

#endif