// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef LIM_SOLVER_2D_H
#define LIM_SOLVER_2D_H

#include "LIMSolver.h"

class TriangleMesh;

class LIMSolver2D : public LIMSolver
{
public:

  virtual ~LIMSolver2D();

  void Init(DeformableMesh* object);
  void Init(TriangleMesh* object);

protected:
  
  TriangleMesh* mesh;
  
  Eigen::Matrix<int,6,Eigen::Dynamic> TriangleVertexIdx;
  Eigen::Matrix<double,Eigen::Dynamic,1> initialNodes;
  
  Eigen::Matrix<double*,1,Eigen::Dynamic> problemHessianCoeffs;
  Eigen::Matrix<double*,6,Eigen::Dynamic> nonFlipHessianCoeffs;
  Eigen::Matrix<double*,21,Eigen::Dynamic> denseHessianCoeffs;
  Eigen::Matrix<double*,1,Eigen::Dynamic> posConstraintsHessianCoeffs;

  LIMSolver2D();

  // abstract functions
  virtual void debugOutput(std::stringstream& info) = 0;
  virtual void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx) = 0;
  virtual double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x) = 0;
  virtual void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad) = 0;
  virtual void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess) = 0;

private:

  // implementation of abstract functions of LIMSolver
  void computeRestPoseFunctionParameters();

  // implementation of abstract functions of NMSolver
  virtual void prepareNMProblemData();
  virtual double computeNMFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
  virtual void computeNMGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad);
  virtual void computeNMHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x);
};

const int NonFlipHessian2DIdx[6][2] = {{1,2},{0,3},{1,4},{3,4},{0,5},{2,5}};

#endif

