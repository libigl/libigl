// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "UniformLaplacian_LIMSolver3D.h"
#include "TetrahedronMesh.h"

#define IGL_HEADER_ONLY
#include "igl/adjacency_matrix.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"

UniformLaplacian_LIMSolver3D::UniformLaplacian_LIMSolver3D()
{
}

UniformLaplacian_LIMSolver3D::~UniformLaplacian_LIMSolver3D()
{
}

void UniformLaplacian_LIMSolver3D::debugOutput(std::stringstream& info)
{	
  std::cout << "LP:" << info.str() << "\n"; 
}

void UniformLaplacian_LIMSolver3D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();
  const int numTets = mesh->Tetrahedra->rows();

  // find connectivity of provided tet mesh
  std::vector<std::pair<int,int> > edges;
  const int tetEdges[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
  
  for(int n=0;n<numTets;n++)
  {
    for(int e=0;e<6;e++)
    {
      int node0 = mesh->Tetrahedra->coeff(n,tetEdges[e][0]);
      int node1 = mesh->Tetrahedra->coeff(n,tetEdges[e][1]);

      if(node0 < node1)
        edges.push_back(std::pair<int,int>(node0,node1));
      else
        edges.push_back(std::pair<int,int>(node1,node0));
    }
  }
  
  // std::sort edges
  std::sort(edges.begin(),edges.end());
  
  // remove dublicates
  std::vector<std::pair<int,int> >::iterator end;
  end = std::unique(edges.begin(),edges.end());

  // create sparse uniform laplacian matrix L
  L.resize(numVariables,numVariables);
  std::vector<Eigen::Triplet<double> > triplets;
  triplets.reserve(numVariables);
  for(std::vector<std::pair<int,int> >::iterator iter=edges.begin();iter!=end;++iter)
  {
    int node0 = iter->first;
    int node1 = iter->second;
    
    for(int i=0;i<3;i++)
    {
      int row = node0*3+i;
      int col = node1*3+i;
      
      triplets.push_back(Eigen::Triplet<double>(col,row,-1));
      triplets.push_back(Eigen::Triplet<double>(row,col,-1));
    }
    
    for(int i=0;i<3;i++)
    {
      int index = node0*3+i;
      triplets.push_back(Eigen::Triplet<double>(index,index,1));
      index = node1*3+i;
      triplets.push_back(Eigen::Triplet<double>(index,index,1));
    }
  }
  L.setFromTriplets(triplets.begin(),triplets.end());

  // bi-harmonic laplacian
  L = L*L;

  TetrahedronVertexIdx.resize(12,mesh->Tetrahedra->rows());

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

double UniformLaplacian_LIMSolver3D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // laplacian energy function f(v) = (v-p)'L(v-p)
  Eigen::VectorXd vv0 = x-initialNodes;
  return Divider * 0.5 * vv0.transpose() * L * vv0;
}

void UniformLaplacian_LIMSolver3D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // laplacian
  grad = L * (x-initialNodes) * Divider;
}

void UniformLaplacian_LIMSolver3D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // laplacian
  int numElem = 0;
  for (int k=0;k<L.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L,k);it;++it)
    {
      if(it.row() <= it.col())
        *hess[numElem++] = it.value() * Divider;
    }
  }
}
