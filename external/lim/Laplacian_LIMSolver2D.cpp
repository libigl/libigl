// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "Laplacian_LIMSolver2D.h"
#include "TriangleMesh.h"

#define IGL_HEADER_ONLY
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"

Laplacian_LIMSolver2D::Laplacian_LIMSolver2D()
{
}

Laplacian_LIMSolver2D::~Laplacian_LIMSolver2D()
{
}

void Laplacian_LIMSolver2D::debugOutput(std::stringstream& info)
{
  std::cout << "LP: " << info.str() << "\n"; 
}

void Laplacian_LIMSolver2D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();

  // create sparse laplacian matrix L
  Eigen::SparseMatrix<double> tempL;
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(*mesh->InitalVertices,*mesh->Triangles,igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::cotmatrix(*mesh->InitalVertices, *mesh->Triangles,tempL);
  tempL *= -1;
  
  // compute inverse of mass matrix
  for (int k=0;k<M.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(M,k);it;++it)
    {
      it.valueRef() = 1.0/it.value();
    }
  }

  tempL = tempL*M*tempL;

  // Create L matrix
  L.resize(numVariables,numVariables);
  std::vector<Eigen::Triplet<double> > triplets;
  for (int k=0;k<tempL.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(tempL,k);it;++it)
    {
      int row = 2*it.row();
      int col = 2*it.col();
      triplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
      triplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
    }
  }
  L.setFromTriplets(triplets.begin(), triplets.end());
  
  for (int k=0;k<L.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,k);it;++it)
    {
      int row = it.row();
      int col = it.col();
      
      // std::sort for upper triangule matrix
      if(row <= col)
      {
        hessRowIdx.push_back(row);
        hessColIdx.push_back(col);
      }
    }
  }
}

double Laplacian_LIMSolver2D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // laplacian energy function f(v) = (v-p)'L(v-p)
  Eigen::VectorXd vv0 = x-initialNodes;
  return 0.5 * vv0.transpose() * L * vv0;
}

void Laplacian_LIMSolver2D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // laplacian
  grad = L * (x-initialNodes);
}

void Laplacian_LIMSolver2D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // laplacian
  int numElem = 0;
  for (int k=0;k<L.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,k);it;++it)
    {
      if(it.row() <= it.col())
        *hess[numElem++] = it.value();
    }
  }
}
