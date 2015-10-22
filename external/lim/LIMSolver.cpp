// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LIMSolver.h"
#include "DeformableMesh.h"

#include <assert.h>

LIMSolver::LIMSolver()
{
  // switches
  EnableBarriers = true;
  EnableNeoHookeanBarriers = false;
  EnableLogBarriers = false;
  EnableBarrierCompensation = false;
  EnableSubstepping = true;
  EnableAlpaUpdate= true;
  EnableOutput = true;

  numIterations= 0;

  // parameters
  AlphaRatio = 1e3;
  Alpha = 1e8;
  Beta = 0.01;
  Gamma = 1;
  CompensationExp = 1;
  Divider = 1;

  // Substepping
  MaxSubStep = 1;
  SubStepExp = 2;

  // output
  CurrentPositionalEnergy = 0;
  CurrentPositionalSubStepEnergy = 0;
  CurrentConstraintEnergy = 0;
}

LIMSolver::~LIMSolver()
{
}

void LIMSolver::Init(DeformableMesh* mesh)
{
  dmesh = mesh;
  
  // check initialization of mesh
  assert(mesh->InitalVertices != NULL);
  assert(mesh->DeformedVertices != NULL);
  assert(mesh->PredictedVertices != NULL);
  //assert(mesh->BorderVertices != NULL);
  assert(mesh->ConstraintMatrix != NULL);
  assert(mesh->ConstraintTargets != NULL);

  if(EnableOutput)
    std::cout << "Initializing energy...";

  UpdatePositionalConstraintMatrix();

  NMSolver::init();

  computeRestPoseFunctionParameters();

  if(EnableOutput)
    std::cout << " done\n";
  
  // init map to positional constraint vector
  //new (&positionalConstraints) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>,Eigen::Aligned,Eigen::Stride<1,3> >(dmesh->PositionalConstraints->data(), dmesh->PositionalConstraints->rows()*dim,1);
}

void LIMSolver::UpdatePositionalConstraintMatrix()
{
  subStepConstraints.resize(dmesh->ConstraintTargets->rows());
  // initialize positional constraints structure
  linearConstraintsMatrix2 = dmesh->ConstraintMatrix->transpose()* *dmesh->ConstraintMatrix;
}

void LIMSolver::debugOutput()
{
  std::stringstream output (std::stringstream::in | std::stringstream::out);
  std::stringstream tempOutput (std::stringstream::in | std::stringstream::out);

  tempOutput << numIterations;
  output << std::left << std::setw(6) << tempOutput.str(); tempOutput.str("");
  tempOutput << "H+" << std::setprecision(3) << CurrentLambda;
  output << std::left << std::setw(10) << tempOutput.str(); tempOutput.str("");
  tempOutput << "LS:" << std::setprecision(3) << CurrentStepSize;
  output << std::left << std::setw(16) << tempOutput.str(); tempOutput.str("");
  tempOutput << "P:" << std::setprecision(6) << CurrentPositionalEnergy;
  output << std::left << std::setw(16) << tempOutput.str(); tempOutput.str("");
  tempOutput << "C:" << std::setprecision(6) << CurrentConstraintEnergy;
  output << std::left << std::setw(16) << tempOutput.str(); tempOutput.str("");
  tempOutput << "E:" << std::setprecision(8) << CurrentDeformationEnergy;
  output << std::left << std::setw(18) << tempOutput.str(); tempOutput.str("");

  debugOutput(output);
}

void LIMSolver::beforeSolve()
{
  // initial solution
  for(int i=0;i<numVariables;i++)
    solution[i] = dmesh->DeformedVertices->coeff(VertexPositionIndices[i]);

  // update substeps
  if(EnableSubstepping == false)
    subStepConstraints = *dmesh->ConstraintTargets;
}

void LIMSolver::afterSolve()
{ 
  // update Alpha
  if(EnableAlpaUpdate)
  {
    double newAlpha = computeAlpha();
    if(newAlpha > Alpha)
      Alpha = newAlpha;
  }
  else
  {
     // compute squared distance of positional constraints
     CurrentPositionalEnergy = (*dmesh->ConstraintMatrix * solution - *dmesh->ConstraintTargets).squaredNorm();
  }

  // apply solution
  for(int i=0;i<numVariables;i++)
  {
    dmesh->PredictedVertices->coeffRef(VertexPositionIndices[i]) = solution[i];
  }
}

double LIMSolver::computeAlpha()
{
  double alpha = 1e16;

  // compute squared distance of positional constraints
  CurrentPositionalEnergy = (*dmesh->ConstraintMatrix * solution - *dmesh->ConstraintTargets).squaredNorm();
  
  if(CurrentPositionalEnergy > 0)
  {
    alpha = std::abs(CurrentDeformationEnergy + CurrentConstraintEnergy) * AlphaRatio / CurrentPositionalEnergy;
    if(alpha < AlphaRatio) alpha = AlphaRatio;
    if(alpha > 1e16) alpha = 1e16;
  }

  return alpha;
}

void LIMSolver::updateSubStepping()
{
  if(CurrentLambda == 0)
  {
    subStepConstraints = *dmesh->ConstraintTargets;
  }
  else
  {
    double subStepSize = (1/(1+pow(CurrentLambda,SubStepExp)));
    subStepConstraints = subStepSize * *dmesh->ConstraintTargets - (subStepSize-1) * *dmesh->ConstraintMatrix * solution;
  }
}

void LIMSolver::afterHessianFactorization()
{
  if(EnableSubstepping)
  {
    this->updateSubStepping();
    functionValue = computeNMFunction(solution);
    computeNMGradient(solution, gradient);
  }
}

void LIMSolver::Restart()
{
  stepSize = 1;
  numIterations = 0;
  if(EnableAlpaUpdate)
  {
    Alpha = computeAlpha();
  }
  else
  {
    // compute squared distance of positional constraints
    CurrentPositionalEnergy = (*dmesh->ConstraintMatrix * solution - *dmesh->ConstraintTargets).squaredNorm();
  }
}

int LIMSolver::Solve()
{
  beforeSolve();
  
  NMSolver::solve();

  afterSolve();

  if(EnableOutput)
    LIMSolver::debugOutput();

  numIterations++;

  int result = 1;
  if(CurrentFV == std::numeric_limits<double>::infinity())
    result = -1;

  return result;
}

void LIMSolver::getNMProblemSize()
{
  numVariables = dmesh->InitalVertices->rows()*dim;
}
