// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef LIM_SOLVER_3D_H
#define LIM_SOLVER_3D_H

#include "LIMSolver.h"

class DeformableMesh;
class TetrahedronMesh;

class LIMSolver3D : public LIMSolver
{
public:

  virtual ~LIMSolver3D();

  void Init(DeformableMesh* object);
  void Init(TetrahedronMesh* object);

protected:
  
  TetrahedronMesh* mesh;

  Eigen::Matrix<int,12,Eigen::Dynamic> TetrahedronVertexIdx;
  Eigen::Matrix<double,Eigen::Dynamic,1> initialNodes;
  
  Eigen::Matrix<double*,1,Eigen::Dynamic> problemHessianCoeffs;
  Eigen::Matrix<double*,36,Eigen::Dynamic> nonFlipHessianCoeffs;
  Eigen::Matrix<double*,78,Eigen::Dynamic> denseHessianCoeffs;
  Eigen::Matrix<double*,1,Eigen::Dynamic> posConstraintsHessianCoeffs;
    
  LIMSolver3D();

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

const int NonFlipHessian3DIdx[36][2] = {{1,3},{2,3},{0,4},{2,4},{0,5},{1,5},{1,6},{2,6},{4,6},{5,6},
                      {0,7},{2,7},{3,7},{5,7},{0,8},{1,8},{3,8},{4,8},{1,9},{2,9},{4,9},{5,9},{7,9},{8,9},
                    {0,10},{2,10},{3,10},{5,10},{6,10},{8,10},{0,11},{1,11},{3,11},{4,11},{6,11},{7,11}};

#endif

