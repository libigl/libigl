// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LGARAP_LIMSolver3D.h"
#include "TetrahedronMesh.h"

#include "igl/svd3x3.h"

#define IGL_HEADER_ONLY
#include "igl/cotmatrix_entries.h"

LGARAP_LIMSolver3D::LGARAP_LIMSolver3D()
{
}

LGARAP_LIMSolver3D::~LGARAP_LIMSolver3D()
{
}

int LGARAP_LIMSolver3D::Solve()
{
  computeLocalStep();

  return LIMSolver3D::Solve();
}

void LGARAP_LIMSolver3D::computeLocalStep()
{
  const int numVertices = mesh->InitalVertices->rows();
  const int numTets = mesh->Tetrahedra->rows();
  
  Eigen::MatrixXd uu, CovMat;
  Eigen::Matrix3d rot;
  Eigen::Matrix<double,6,6> cc;
  cc.fill(0);
  uu.resize(6,3);
  CovMat.resize(3,3);
  Eigen::Vector3d u[4];

  R.resize(9*numTets);
      
  // local step: Compute best rigid transformations
  for (int t=0;t<numTets;t++)
  {
    Eigen::Matrix<int,4,1> indices = mesh->Tetrahedra->row(t);
      
    for(int i=0;i<4;i++)
      u[i] = mesh->DeformedVertices->row(indices[i]);

    for(int i=0;i<6;i++)
      uu.row(i) = u[TetEdgeVertices[i][1]] - u[TetEdgeVertices[i][0]];

    for(int i=0;i<6;i++)
      cc(i,i) = CotanWeights(t,i);
    
    CovMat = RestPoseEdges.block<3,6>(3*t,0) * cc * uu;

    Eigen::Matrix3f A = CovMat.cast<float>();
    Eigen::Matrix<float,3,3> U, Vt;
    Eigen::Matrix<float,3,1> S;
    
    igl::svd3x3(A, U, S, Vt);
    rot = (Vt * U.transpose()).cast<double>();

    const int idx = 9*t;
    for(int x=0;x<3;x++)
      for(int y=0;y<3;y++)
        R(idx+x*3+y) = rot(y,x);
  }
}

void LGARAP_LIMSolver3D::debugOutput(std::stringstream& info)
{
  std::cout << "AR: " << info.str() << "\n"; 
}

void LGARAP_LIMSolver3D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numVertices = mesh->InitalVertices->rows();
  const int numTets = mesh->Tetrahedra->rows();

  CotanWeights.resize(numTets,6);
  igl::cotmatrix_entries(*mesh->InitalVertices,*mesh->Tetrahedra,CotanWeights);

  // Create matrices L, K
  Eigen::SparseMatrix<double> B, tempL, tempK, restV;
  tempL.resize(numVertices,numVertices);
  tempK.resize(3*numTets,numVertices);
  std::vector<Eigen::Triplet<double> > LTriplets, BTriplets, KTriplets, VTriplets;
  for(int t=0;t<numTets;t++)
  {
    B.resize(numVertices,6);
    restV.resize(numVertices,3);
    BTriplets.clear();
    VTriplets.clear();
    
    Eigen::Matrix<int,4,1> indices = mesh->Tetrahedra->row(t);
    for(int i=0;i<6;i++)
    {
      int vIdx0 = indices(TetEdgeVertices[i][0]);
      int vIdx1 = indices(TetEdgeVertices[i][1]);

      // Create incident matrix B_i for i'th triangle
      BTriplets.push_back(Eigen::Triplet<double>(vIdx0,i,1));
      BTriplets.push_back(Eigen::Triplet<double>(vIdx1,i,-1));

      // Create 3D tet rest pose vertex matrix
      Eigen::Vector3d v0 = mesh->InitalVertices->row(vIdx0);
      Eigen::Vector3d v1 = mesh->InitalVertices->row(vIdx1);
      VTriplets.push_back(Eigen::Triplet<double>(vIdx0,0,v0(0)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx0,1,v0(1)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx0,2,v0(2)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx1,0,v1(0)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx1,1,v1(1)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx1,2,v1(2)));

      // cotangent gives not cot but multiple of opposite edge length l: l/6*cot
      CotanWeights(t,i) = CotanWeights(t,i)*6/(v1-v0).norm();
    }
    B.setFromTriplets(BTriplets.begin(),BTriplets.end());
    restV.setFromTriplets(VTriplets.begin(),VTriplets.end());

    // cotangent weights
    Eigen::SparseMatrix<double> Cm(6,6);
    for(int i=0;i<6;i++)
      Cm.insert(i,i) = CotanWeights(t,i);

    // Create B*C*B'
    Eigen::SparseMatrix<double> BCBT = B*Cm*B.transpose();

     // Stack up K temp matrix
    Eigen::SparseMatrix<double> tempK = restV.transpose()*BCBT/3.0;
    for (int k=0;k<tempK.outerSize();++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(tempK,k);it;++it)
      {
        KTriplets.push_back(Eigen::Triplet<double>(3*t+it.row(),it.col(),it.value()));
      }
    }

    // Sum up L temp matrix
    for (int k=0;k<BCBT.outerSize();++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(BCBT,k);it;++it)
      {
        int row = 3*it.row();
        int col = 3*it.col();
        LTriplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
        LTriplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
        LTriplets.push_back(Eigen::Triplet<double>(row+2,col+2,it.value()));
      }
    }
  }

   // Create L matrix
  L.resize(numVariables,numVariables);
  L.setFromTriplets(LTriplets.begin(), LTriplets.end());

  // Create K matrix
  tempK.setFromTriplets(KTriplets.begin(),KTriplets.end());
  K.resize(tempK.rows()*3,tempK.cols()*3);
  KTriplets.clear();
  for (int k=0;k<tempK.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(tempK,k);it;++it)
    {
      int row = 3*it.row();
      int col = 3*it.col();
      KTriplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
      KTriplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
      KTriplets.push_back(Eigen::Triplet<double>(row+2,col+2,it.value()));
    }
  }
  K.setFromTriplets(KTriplets.begin(), KTriplets.end());

  RestPoseEdges.resize(3*numTets,6);
  Eigen::Vector3d p[4];
  for (int t=0;t<numTets;t++)
  {
    Eigen::Matrix<int,4,1> indices = mesh->Tetrahedra->row(t);
    
    for(int i=0;i<4;i++)
      p[i] = mesh->InitalVertices->row(indices[i]);

    for(int i=0;i<6;i++)
      RestPoseEdges.block<3,1>(3*t,i) = p[TetEdgeVertices[i][1]] - p[TetEdgeVertices[i][0]];
  }

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

  // compute constant energy bias for restpose
  int numVariables = numVertices*3;
  Eigen::Matrix<double,Eigen::Dynamic,1> restPose(numVariables);
  for(int n=0;n<numVertices;n++)
  {
    for(int i=0;i<3;i++)
      restPose[n*3+i] = mesh->InitalVertices->coeff(n,i);
  }
  constantEnergyPart = 0.5*restPose.transpose()*L*restPose;
}

double LGARAP_LIMSolver3D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // ARAP global step energy
  double xTLx = x.transpose()*L*x;
  double KTRx = R.transpose()*K*x;
  return 0.5*xTLx - KTRx + constantEnergyPart;
}

void LGARAP_LIMSolver3D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // ARAP global step
  grad = L*x - K.transpose()*R;
}

void LGARAP_LIMSolver3D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // ARAP global step
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
