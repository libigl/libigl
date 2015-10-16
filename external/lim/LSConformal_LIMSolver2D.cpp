// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LSConformal_LIMSolver2D.h"
#include "TriangleMesh.h"

#define IGL_HEADER_ONLY
#include "igl/cotmatrix.h"

LSConformal_LIMSolver2D::LSConformal_LIMSolver2D()
{
}

LSConformal_LIMSolver2D::~LSConformal_LIMSolver2D()
{
}

void LSConformal_LIMSolver2D::debugOutput(std::stringstream& info)
{
  std::cout << "LP: " << info.str() << "\n"; 
}

void LSConformal_LIMSolver2D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();
  const int numTriangles = mesh->Triangles->rows();

  // create sparse laplacian matrix L
  Eigen::SparseMatrix<double> tempL;//(numNodes,numNodes);
  igl::cotmatrix(*mesh->InitalVertices, *mesh->Triangles,tempL);
  tempL *= -1;

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

  // Create A matrix
  int numBorderVertices = mesh->BorderVertices->rows();
  
  Eigen::SparseMatrix<double> B, C, BNonZero, CNonZero;
  A.resize(2*numNodes,2*numNodes);
  B.resize(2*numBorderVertices, 2*numNodes);
  C.resize(2*numBorderVertices, 2*numNodes);
  
  std::vector<Eigen::Triplet<double> > tripletsB, tripletsC;
  for(int i=0,j=numBorderVertices-1;i<mesh->BorderVertices->rows();j=i,i++)
  {
    int rIdx = 2*i;
    int iIdx = 2*mesh->BorderVertices->coeff(i); 	
    int jIdx = 2*mesh->BorderVertices->coeff(j);

    tripletsB.push_back(Eigen::Triplet<double>(rIdx,jIdx,1));
    tripletsB.push_back(Eigen::Triplet<double>(rIdx+1,jIdx+1,-1));
    tripletsC.push_back(Eigen::Triplet<double>(rIdx,iIdx+1,1));
    tripletsC.push_back(Eigen::Triplet<double>(rIdx+1,iIdx,1));
  }
  B.setFromTriplets(tripletsB.begin(), tripletsB.end());
  C.setFromTriplets(tripletsC.begin(), tripletsC.end());

  A = B.transpose()*C;
  Eigen::SparseMatrix<double> AT = A.transpose();
  L = L - 0.5*(A+AT);

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

double LSConformal_LIMSolver2D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // laplacian energy function f(v) = (v-p)'L(v-p)
  return 0.5 * x.transpose() * L * x;
}

void LSConformal_LIMSolver2D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // laplacian
  grad = L * x;
}

void LSConformal_LIMSolver2D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
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
