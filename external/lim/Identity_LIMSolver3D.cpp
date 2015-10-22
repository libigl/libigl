// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "Identity_LIMSolver3D.h"
#include "TetrahedronMesh.h"

Identity_LIMSolver3D::Identity_LIMSolver3D()
{
}

Identity_LIMSolver3D::~Identity_LIMSolver3D()
{
}

void Identity_LIMSolver3D::debugOutput(std::stringstream& info)
{
  std::cout << "I: " << info.str() << "\n"; 
}

void Identity_LIMSolver3D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();
  const int numTets = mesh->Tetrahedra->rows();

  // create sparse identity matrix
  Identity.resize(numVariables,numVariables);
  std::vector<Eigen::Triplet<double> > triplets;
  triplets.reserve(numVariables);
  for(int i=0;i<numVariables;i++)
  {
    triplets.push_back(Eigen::Triplet<double>(i,i,1));	
  }
  Identity.setFromTriplets(triplets.begin(),triplets.end());

  for (int i=0;i<numVariables;i++)
  {
    hessRowIdx.push_back(i);
    hessColIdx.push_back(i);
  }
}

double Identity_LIMSolver3D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // laplacian energy function f(v) = (v-p)'L(v-p)
  Eigen::VectorXd vv0 = x-initialNodes;
  return 0.5 * vv0.transpose() * Identity * vv0;
}

void Identity_LIMSolver3D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // laplacian
  grad = x-initialNodes;
}

void Identity_LIMSolver3D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // identity
  for(int i=0;i<numVariables;i++)
    *hess[i] = 1;
}
