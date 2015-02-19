// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef LG_ARAP_LIM_SOLVER_2D_H
#define LG_ARAP_LIM_SOLVER_2D_H

#include "LIMSolver2D.h"

class LGARAP_LIMSolver2D : public LIMSolver2D
{
public:

  LGARAP_LIMSolver2D();
  ~LGARAP_LIMSolver2D();

  int Solve();

private:
  
  Eigen::SparseMatrix<double> L, K;
  Eigen::Matrix<double,Eigen::Dynamic,1> R;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> CotanWeights;
  Eigen::Matrix<double,Eigen::Dynamic,3> RestPoseEdges;
  double constantEnergyPart;

  void computeLocalStep();

  // implementation of abstract functions of LIMSolver2D
  void debugOutput(std::stringstream& info);
  void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx);
  double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess);
};

// edge order in igl::cotangent
const int TriEdgeVertices[3][2] = {{1,2},{2,0},{0,1}};

#endif