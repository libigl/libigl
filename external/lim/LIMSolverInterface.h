// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

//----------------------------------------------------------------------------------------
// LIMSolverInterface.h
// Date: 07.06.13
// Author: Christian Schüller
//----------------------------------------------------------------------------------------

#pragma once

#ifndef LIM_SOLVER_INTERFACE_H
#define LIM_SOLVER_INTERFACE_H

#include "TriangleMesh.h"
#include "TetrahedronMesh.h"

#include "LIMSolver2D.h"
#include "LIMSolver3D.h"

#include "Identity_LIMSolver2D.h"
#include "Dirichlet_LIMSolver2D.h"
#include "UniformLaplacian_LIMSolver2D.h"
#include "Laplacian_LIMSolver2D.h"
#include "GreenStrain_LIMSolver2D.h"
#include "LGARAP_LIMSolver2D.h"
#include "LSConformal_LIMSolver2D.h"
#include "Poisson_LIMSolver2D.h"

#include "Identity_LIMSolver3D.h"
#include "Dirichlet_LIMSolver3D.h"
#include "UniformLaplacian_LIMSolver3D.h"
#include "Laplacian_LIMSolver3D.h"
#include "GreenStrain_LIMSolver3D.h"
#include "LGARAP_LIMSolver3D.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

// LIM data structure
struct LIMData
{
  bool isTetMesh;
  DeformableMesh* mesh;
  LIMSolver* solver;
  int iteration;
};

//----------------------------------------------------------------------------------------
// Function: FreeLIMData
//----------------------------------------------------------------------------------------
// Description:
// Releases given LIM data object
//----------------------------------------------------------------------------------------
// Input:
// data             Pointer to LIMData instance
//----------------------------------------------------------------------------------------
void FreeLIMData(LIMData* data)
{
  delete data->mesh->InitalVertices;
  delete data->mesh->DeformedVertices;
  delete data->mesh->PredictedVertices;
  
  if(data->isTetMesh)
    delete static_cast<TetrahedronMesh*>(data->mesh)->Tetrahedra;
  else
    delete static_cast<TriangleMesh*>(data->mesh)->Triangles;

  delete data->mesh->BorderVertices;
  delete data->mesh->ConstraintMatrix;
  delete data->mesh->ConstraintTargets;
  
  delete data->mesh;
  delete data->solver;
}

//----------------------------------------------------------------------------------------
// Function: InitLIM
//----------------------------------------------------------------------------------------
// Description:
// Initializes the LIM Solver data before calling the function ComputeLIM_Step
//----------------------------------------------------------------------------------------
// Input:
// vertices          vx3 matrix containing vertex position of the mesh
// initialVertices   vx3 matrix containing vertex position of initial rest pose mesh
// elements          exd matrix containing vertex indices of all elements
// borderVertices    (optional) only needed for 2D LSCM) vector containing indices of border vertices
// gradients         (optional) only needed for 2D Poisson) vector containing partial derivatives of target element gradients (structure is: [xx_1, xy_1, xx_2, xy_2, ..., xx_v, xy_v, yx_1, yy_1, yx_2, yy_2, ..., yx_v, yy_v]')
// constraintMatrix  C: (c)x(3xv) sparse linear positional constraint matrix
//                   X,Y,Z-coordinates are alternatingly stacked per row (structure for triangles: [x_1, y_1, z_1, x_2, y_2, z_2, ..., x_v,y_v,z_v])
//                   and each row of C belongs to a linear constraint.
// constraintTargets d: c vector target positions
// energyType        type of used energy: 0=Dirichlet,1=Laplacian,2=Green,3=ARAP,4=LSCM,5=Poisson
// enableOutput      (optional) enables the output (#iteration / hessian correction / step size / positional constraints squared error / barrier constraints energy / deformation energy)
// enableBarriers    (optional) enables the non-flip constraints (default = true)
// enableAlphaUpdate (optional) enables dynamic alpha weight adjustment (default = true)
// beta              (optional) steepness factor of barrier slopes (default: ARAP/LSCM = 0.01, Green = 1)
// eps               (optional) smallest valid triangle area (default: 1e-5 * smallest triangle)
//
// where:
// v : # vertices
// c : # linear constraints
// e : # elements of mesh
// d : # vetices per element (triangle = 3, tet = 4)
//----------------------------------------------------------------------------------------
// Return value:
// data              a pointer to the LIM data object
//----------------------------------------------------------------------------------------
LIMData* InitLIM(
  Eigen::Matrix<double,Eigen::Dynamic,3>& vertices,
  const Eigen::Matrix<double,Eigen::Dynamic,3>& initialVertices,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>& elements,
  const std::vector<int>& borderVertices,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& gradients,
  const Eigen::SparseMatrix<double>& constraintMatrix,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& constraintTargets,
  int energyType,
  bool enableOuput = true,
  bool enableBarriers = true,
  bool enableAlphaUpdate = true,
  double beta = -1,
  double eps = -1)
{
  LIMData* data = new LIMData();
  data->isTetMesh = (elements.cols() == 4);

  //------------------------------------------------------------------------------------
  // Init mesh object
  //------------------------------------------------------------------------------------
  DeformableMesh* mesh = NULL;

  if(data->isTetMesh)
  {
    TetrahedronMesh* tetMesh = new TetrahedronMesh();
    mesh = tetMesh;

    tetMesh->Tetrahedra = new Eigen::Matrix<int,Eigen::Dynamic,4>(elements);
  }
  else
  {
    TriangleMesh* triMesh = new TriangleMesh();
    mesh = triMesh;

    triMesh->Triangles = new Eigen::Matrix<int,Eigen::Dynamic,3>(elements);
    triMesh->BorderVertices = new Eigen::Matrix<int,Eigen::Dynamic,1>();
    triMesh->IsCorotatedTriangles = false;

    triMesh->BorderVertices->resize(borderVertices.size(),1);
    for(int i=0;i<(int)borderVertices.size();i++)
      triMesh->BorderVertices->coeffRef(i) = borderVertices[i];
  }
  
  mesh->InitalVertices = new Eigen::Matrix<double,Eigen::Dynamic,3>(initialVertices);
  mesh->DeformedVertices = new Eigen::Matrix<double,Eigen::Dynamic,3>(vertices);
  mesh->PredictedVertices = new Eigen::Matrix<double,Eigen::Dynamic,3>(vertices);  
  mesh->ConstraintMatrix = new Eigen::SparseMatrix<double>(constraintMatrix);
  mesh->ConstraintTargets = new Eigen::Matrix<double,Eigen::Dynamic,1>(constraintTargets);

  mesh->InitMesh();

  if(eps != -1) mesh->EPS3 = eps;
  
  //------------------------------------------------------------------------------------
  // Intit solver
  //------------------------------------------------------------------------------------
  
  LIMSolver* solver = NULL;
  
  if(data->isTetMesh)
  {
    switch(energyType)
    {
      case 0:
        solver = new Dirichlet_LIMSolver3D();
        break;

      case 1:
        solver = new Laplacian_LIMSolver3D();
        break;

      case 2:
        solver = new GreenStrain_LIMSolver3D();
        break;

      case 3:
        solver = new LGARAP_LIMSolver3D();
        break;
  
      default:
        solver = new GreenStrain_LIMSolver2D();
        break;
    }
  }
  else
  {
    switch(energyType)
    {
      case 0:
        solver = new Dirichlet_LIMSolver2D();
        break;

      case 1:
        solver = new Laplacian_LIMSolver2D();
        break;

      case 2:
        solver = new GreenStrain_LIMSolver2D();
        break;

      case 3:
        solver = new LGARAP_LIMSolver2D();
        break;

      case 4:
        solver = new LSConformal_LIMSolver2D();
        break;

      case 5:
      {
        Poisson_LIMSolver2D* psolver = new Poisson_LIMSolver2D();
        psolver->b = gradients;
        solver = psolver;
      }
      break;
  
      default:
        new GreenStrain_LIMSolver2D();
        break;
    }
  }
  
  solver->Init(mesh);

  solver->EnableBarriers = enableBarriers;
  if(beta != -1) solver->Beta = beta;

  data->mesh = mesh;
  data->solver = solver;
  data->iteration = 0;

  return data;
}

LIMData* InitLIM(
  Eigen::Matrix<double,Eigen::Dynamic,3>& vertices,
  const Eigen::Matrix<double,Eigen::Dynamic,3>& initialVertices,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>& elements,
  const Eigen::SparseMatrix<double>& constraintMatrix,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& constraintTargets,
  int energyType,
  bool enableOuput = true,
  bool enableBarriers = true,
  bool enableAlphaUpdate = true,
  double beta = -1,
  double eps = -1)
{
  vector<int> borderVertices;
  Eigen::VectorXd gradients;

  return InitLIM(
    vertices,
    initialVertices,
    elements,
    borderVertices,
    gradients,
    constraintMatrix,
    constraintTargets,
    energyType,
    enableOuput,
    enableBarriers,
    enableAlphaUpdate,
    beta,
    eps);
}

//----------------------------------------------------------------------------------------
// Function: ComputeLIM
//----------------------------------------------------------------------------------------
// Description:
// Computes a locally injective mapping of a triangle or tet-mesh based on a deformation energy
// subject to some provided linear positional constraints Cv-d.
//----------------------------------------------------------------------------------------
// Input:
// vertices          vx3 matrix containing vertex position of the mesh
// initialVertices   vx3 matrix containing vertex position of initial rest pose mesh
// elements          exd matrix containing vertex indices of all elements
// borderVertices    (optional) (only needed for 2D LSCM) vector containing indices of border vertices
// gradients         (optional) (only needed for 2D Poisson) vector containing partial derivatives of target element gradients (structure is: [xx_1, xy_1, xx_2, xy_2, ..., xx_v, xy_v, yx_1, yy_1, yx_2, yy_2, ..., yx_v, yy_v]')
// constraintMatrix  C: (c)x(3xv) sparse linear positional constraint matrix
//                   X,Y,Z-coordinates are alternatingly stacked per row (structure for triangles: [x_1, y_1, z_1, x_2, y_2, z_2, ..., x_v,y_v,z_v])
//                   and each row of C belongs to a linear constraint.
// constraintTargets d: c vector target positions
// energyType        type of used energy: 0=Dirichlet,1=Laplacian,2=Green,3=ARAP,4=LSCM
// tolerance         max squared positional constraints error
// maxIteration      max number of iterations
// findLocalMinima   iterating until a local minima is found. If not enabled only tolerance must be fulfilled.
// enableOutput      (optional) enables the output (#itaration / hessian correction / step size / positional constraints / barrier constraints / deformation energy) (default : true)
// enableBarriers    (optional) enables the non-flip constraints (default = true)
// enableAlphaUpdate (optional) enables dynamic alpha weight adjustment (default = true)
// beta              (optional) steepness factor of barrier slopes (default: ARAP/LSCM = 0.01, Green = 1)
// eps               (optional) smallest valid triangle area (default: 1e-5 * smallest triangle)
//
// where:
// v : # vertices
// c : # linear constraints
// e : # elements of mesh
// d : # vetices per element (triangle = 3, tet = 4)
//----------------------------------------------------------------------------------------
// Output:
// vertices          vx3 matrix containing resulting vertex position of the mesh
//----------------------------------------------------------------------------------------
// Return values:
//  1 : Successful optimization with fulfilled tolerance
// -1 : Max iteration reached before tolerance was fulfilled
// -2 : not feasible -> has inverted elements (may want to decrease eps?)
//----------------------------------------------------------------------------------------
int ComputeLIM(
  Eigen::Matrix<double,Eigen::Dynamic,3>& vertices,
  const Eigen::Matrix<double,Eigen::Dynamic,3>& initialVertices,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>& elements,
  const std::vector<int>& borderVertices,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& gradients,
  const Eigen::SparseMatrix<double>& constraintMatrix,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& constraintTargets,
  int energyType,
  double tolerance,
  int maxIteration,
  bool findLocalMinima,
  bool enableOuput = true,
  bool enableBarriers = true,
  bool enableAlphaUpdate = true,
  double beta = -1,
  double eps = -1)
{
  LIMData* data = InitLIM(vertices, initialVertices, elements, borderVertices, gradients, constraintMatrix, constraintTargets, energyType, enableOuput, enableBarriers, enableAlphaUpdate, beta, eps);

  int result = 0;
  while(result == 0)
  {
    if(data->solver->CurrentStepSize < 1e-15 || (data->solver->CurrentPositionalEnergy <= tolerance && (findLocalMinima == false || data->solver->CurrentStepSize < 1e-15)))
      result = 1; // termination criteria fulfilled

    if(data->iteration >= maxIteration)
      result = -1; // max iteration reached
    
    if(result == 0)
    {
      if(data->solver->Solve() == -1)
        result = -2; // state not feasible -> inverted elements
      else
      {
        // swap vertex buffers
        Eigen::Matrix<double,Eigen::Dynamic,3>* temp = data->mesh->DeformedVertices;
        data->mesh->DeformedVertices = data->mesh->PredictedVertices;
        data->mesh->PredictedVertices = temp;

        data->iteration++;
      }
    }
  }

  // assign resulting vertices
  vertices = *data->mesh->DeformedVertices;

  // release solver data
  FreeLIMData(data);

  return result;
}

int ComputeLIM(
  Eigen::Matrix<double,Eigen::Dynamic,3>& vertices,
  const Eigen::Matrix<double,Eigen::Dynamic,3>& initialVertices,
  const Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>& elements,
  const Eigen::SparseMatrix<double>& constraintMatrix,
  const Eigen::Matrix<double,Eigen::Dynamic,1>& constraintTargets,
  int energyType,
  double tolerance,
  int maxIteration,
  bool findLocalMinima,
  bool enableOuput = true,
  bool enableBarriers = true,
  bool enableAlphaUpdate = true,
  double beta = -1,
  double eps = -1)
{
  vector<int> borderVertices;
  Eigen::VectorXd gradients;

  return ComputeLIM(
    vertices,
    initialVertices,
    elements,
    borderVertices,
    gradients,
    constraintMatrix,
    constraintTargets,
    energyType,
    tolerance,
    maxIteration,
    findLocalMinima,
    enableOuput,
    enableBarriers,
    enableAlphaUpdate,
    beta,
    eps);
}

//----------------------------------------------------------------------------------------
// Function: ComputeLIM_Step
//----------------------------------------------------------------------------------------
// Description:
// Computes one minimization step for the given LIM problem. Use InitLim to initialize LIM data.
//----------------------------------------------------------------------------------------
// Input:
// data              LIM data structure
//----------------------------------------------------------------------------------------
// Output:
// vertices          vx3 matrix containing resulting vertex position of the mesh
//----------------------------------------------------------------------------------------
// Return values:
//  1 : Successful optimization step
// -1 : Lim data is not initialized
// -2 : not feasible -> has inverted elements (may want to decrease eps?)
//----------------------------------------------------------------------------------------
int ComputeLIM_Step(
  LIMData*& data,
  Eigen::Matrix<double,Eigen::Dynamic,3>& vertices)
{
  if(data == NULL)
  {
    cerr << "LIM data is not initialized." << endl;
    return -1;
  }
  
  int result = 0;
    
  if(data->solver->Solve() == -1)
    result = -2; // state not feasible -> inverted elements
  else
  {
    // swap vertex buffers
    Eigen::Matrix<double,Eigen::Dynamic,3>* temp = data->mesh->DeformedVertices;
    data->mesh->DeformedVertices = data->mesh->PredictedVertices;
    data->mesh->PredictedVertices = temp;

    data->iteration++;
  }

  // assign resulting vertices
  vertices = *data->mesh->DeformedVertices;

  return result;
}

#endif
