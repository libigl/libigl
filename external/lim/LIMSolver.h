// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef LIM_SOLVER_H
#define LIM_SOLVER_H

#include "NMSolver.h"

class DeformableMesh;

class LIMSolver : public NMSolver
{
public:

  // switches
  bool EnableBarriers;
  bool EnableNeoHookeanBarriers;
  bool EnableLogBarriers;
  bool EnableBarrierCompensation;
  bool EnableSubstepping;
  bool EnableAlpaUpdate;
  bool EnableOutput;

  int numIterations;
  
  // parameters
  double AlphaRatio;
  double Alpha;
  double Beta;
  double Gamma;
  double CompensationExp;
  double Divider;
  
  // substepping
  double MaxSubStep;
  double SubStepExp;

  // output
  double CurrentPositionalEnergy;
  double CurrentPositionalSubStepEnergy;
  double CurrentConstraintEnergy;
  double CurrentDeformationEnergy;

  virtual ~LIMSolver();

  // Update linear positional constrain matrix (structur has to remain the same!)
  void UpdatePositionalConstraintMatrix();

  virtual void Init(DeformableMesh* object);
  virtual int Solve();
  virtual void Restart();

protected:
  DeformableMesh* dmesh;
  int dim;

  Eigen::SparseMatrix<double> linearConstraintsMatrix2;
  Eigen::Matrix<double,Eigen::Dynamic,1> subStepConstraints;

  Eigen::VectorXi VertexPositionIndices;
  Eigen::Matrix<double,Eigen::Dynamic,1> initalSize;
  Eigen::Matrix<double,Eigen::Dynamic,1> barrierParam1;
  Eigen::Matrix<double,Eigen::Dynamic,1> barrierParam2;

  LIMSolver();

  virtual void beforeSolve();
  virtual void afterSolve();
  virtual void debugOutput();

  void updateSubStepping();
  double computeAlpha();
  void computeAlphas(Eigen::Matrix<double,Eigen::Dynamic,1>& alphas);
  
  // override function of NMSolver
  void afterHessianFactorization();

  // implement abstract functions of NMSolver
  void getNMProblemSize();

  // abstract functions of LIMSolver
  virtual void debugOutput(std::stringstream& info) = 0;
  virtual void computeRestPoseFunctionParameters() = 0;
  
  // abstract functions of NMSolver
  virtual void prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx) = 0;
  virtual double computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x) = 0;
  virtual void computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad) = 0;
  virtual void computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess) = 0;

private:

};

#endif
