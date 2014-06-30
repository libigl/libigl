// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef NM_OPTIMIZER_H
#define NM_OPTIMIZER_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class NMSolver
{
public:

  bool EnableWolfeConditions;
  
  double CurrentLambda;
  double CurrentStepSize;
  double CurrentFV;

protected:
  
  bool isPrevHessianSingular;

  // newton step parameters
  double stepSize;
  double minStepSize;
  double maxStepSize;
  
  // hessian correction parameters
  double lambda;
  double minLambda;
  double maxLambda;

  // Wolfe condition step size parameters
  double wc1, wc2;
  double wAlphaMax;

  int numVariables;
  int numNonZeroHessian;

  double functionValue;

  Eigen::SparseMatrix<double> hessian;
  Eigen::Matrix<double,Eigen::Dynamic,1> gradient;
  Eigen::Matrix<double,Eigen::Dynamic,1> solution;
  Eigen::Matrix<double,Eigen::Dynamic,1> prevSolution;
  Eigen::Matrix<double,Eigen::Dynamic,1> step;

  Eigen::Matrix<double*,Eigen::Dynamic,1>  diagHessianCoeffs;

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> lltSolver;

  NMSolver();

  int init();
  int solve();
  bool subSolve();
  double computeStepSizeWolfe(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, double of_x, const Eigen::Matrix<double,Eigen::Dynamic,1>& dir);
  double nocZoomWolfe(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double,Eigen::Dynamic,1>& dir, double slope0, double alphaLo, double alphaHi, double of_0, double ofLo, double c1, double c2);
  void validateWithFD(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  
  // functions to overwrite
  virtual void afterHessianFactorization();

  // abstract functions
  virtual void getNMProblemSize() = 0;
  virtual void prepareNMProblemData() = 0;
  virtual double computeNMFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x) = 0;
  virtual void computeNMGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad) = 0;
  virtual void computeNMHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x) = 0;
};

#endif