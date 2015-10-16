// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "Poisson_LIMSolver2D.h"
#include "TriangleMesh.h"

#define IGL_HEADER_ONLY
#include "igl/grad.h"
#include "igl/doublearea.h"

Poisson_LIMSolver2D::Poisson_LIMSolver2D()
{
  constantEnergyPart = 0;
}

Poisson_LIMSolver2D::~Poisson_LIMSolver2D()
{
}

void Poisson_LIMSolver2D::debugOutput(std::stringstream& info)
{
  std::cout << "LP: " << info.str() << "\n"; 
}
  
void Poisson_LIMSolver2D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();

  // create sparse gradient operator matrix
  Eigen::SparseMatrix<double> tempG;
  Eigen::VectorXd dAreas,dAreasTemp;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> vertices(*mesh->DeformedVertices);
  Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> faces(*mesh->Triangles);
  
  igl::grad(vertices,faces,tempG);

  // Only get x and y derivatives of elements as z is zero
  int newRowSize = 2.0/3.0*tempG.rows();
  std::vector<Eigen::Triplet<double> > triplets;
  for (int k=0;k<tempG.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(tempG,k);it;++it)
    {
      int row = it.row();
      int col = it.col();
      if(row < newRowSize)
      {
        triplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
      }
    }
  }
  tempG.setZero();
  tempG.resize(newRowSize,tempG.cols());
  tempG.setFromTriplets(triplets.begin(), triplets.end());

  // Extend gradient operator matrix for x and y scalar function
  triplets.clear();
  G.resize(newRowSize*2,tempG.cols()*2);
  for (int k=0;k<tempG.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(tempG,k);it;++it)
    {
      int row = it.row()*2;
      int col = it.col()*2;
      triplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
      triplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
    }
  }
  G.setFromTriplets(triplets.begin(), triplets.end());

  // Compute area weights
  Eigen::SparseMatrix<double> M;
  igl::doublearea(vertices,faces,dAreas);
  triplets.clear();
  M.resize(dAreas.rows()*4,dAreas.rows()*4);
  for(int r=0;r<dAreas.rows();r++)
  {
    int id = 4*r;
    triplets.push_back(Eigen::Triplet<double>(id,id,dAreas(r)));
    triplets.push_back(Eigen::Triplet<double>(id+1,id+1,dAreas(r)));
    triplets.push_back(Eigen::Triplet<double>(id+2,id+2,dAreas(r)));
    triplets.push_back(Eigen::Triplet<double>(id+3,id+3,dAreas(r)));
  }
  M.setFromTriplets(triplets.begin(),triplets.end());

  // Compute laplacian
  L = 0.5*G.transpose()*M*G;

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

  GTb = 0.5*G.transpose()*M*b;
  constantEnergyPart = b.transpose()*b;
}

double Poisson_LIMSolver2D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // poisson energy function f(x) = 0.5*||Gx-b||^2 = 0.5*x'Lx - b'Gx + b'b
  return 0.5 * x.transpose() * L * x - GTb.dot(x) + constantEnergyPart;
}

void Poisson_LIMSolver2D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // poisson gradient
  grad = L * x - GTb;
}

void Poisson_LIMSolver2D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // poisson hessian
  int numElem = 0;
  for (int k=0;k<L.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,k);it;++it)
    {
      int row = it.row();
      int col = it.col();

      if(row <= col)
      {
        *hess[numElem++] = it.value();
      }
    }
  }
}
