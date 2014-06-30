// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef GREEN_STRAIN_LIM_SOLVER_3D_H
#define GREEN_STRAIN_LIM_SOLVER_3D_H

#include "LIMSolver3D.h"

class GreenStrain_LIMSolver3D : public LIMSolver3D
{
public:

  GreenStrain_LIMSolver3D();
  ~GreenStrain_LIMSolver3D();
  
private:
  
  Eigen::Matrix<double,4,Eigen::Dynamic> Ms;
  Eigen::Matrix<double,4,Eigen::Dynamic> MMTs;

  // implementation of abstract functions of LIMSolver2D
  void debugOutput(std::stringstream& info);
  void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx);
  double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess);
};

#endif