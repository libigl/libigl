// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LGARAP_LIMSolver2D.h"
#include "TriangleMesh.h"

#define IGL_HEADER_ONLY
#include "igl/cotmatrix_entries.h"

LGARAP_LIMSolver2D::LGARAP_LIMSolver2D()
{
}

LGARAP_LIMSolver2D::~LGARAP_LIMSolver2D()
{
}

int LGARAP_LIMSolver2D::Solve()
{
  computeLocalStep();

  return LIMSolver2D::Solve();
}

void LGARAP_LIMSolver2D::computeLocalStep()
{
  const int numVertices = mesh->InitalVertices->rows();
  const int numTriangles = mesh->Triangles->rows();
  
  Eigen::MatrixXd uu, CovMat;
  Eigen::Matrix2d rot;
  Eigen::Matrix3d cc;
  cc.fill(0);
  uu.resize(3,2);
  CovMat.resize(2,2);
  Eigen::RowVector2d u0,u1,u2;
  Eigen::Vector2d u[3];

  R.resize(4*numTriangles);
      
  // local step: Compute best rigid transformations
  for (int t=0;t<numTriangles;t++)
  {
    Eigen::Vector3i indices = mesh->Triangles->row(t);

    for(int i=0;i<3;i++)
      u[i] = mesh->DeformedVertices->block<1,2>(indices(i),0);

    for(int i=0;i<3;i++)
      uu.row(i) = u[TriEdgeVertices[i][1]]-u[TriEdgeVertices[i][0]];
    
    cc(0,0) = CotanWeights(t,0);
    cc(1,1) = CotanWeights(t,1);
    cc(2,2) = CotanWeights(t,2);
    
    CovMat = RestPoseEdges.block<2,3>(2*t,0) * cc * uu;

    Eigen::JacobiSVD<Eigen::MatrixXd> svdOfCovMat(CovMat, Eigen::ComputeThinU | Eigen::ComputeThinV);
        
    Eigen::Matrix2d u = svdOfCovMat.matrixU();
    Eigen::Vector2d s = svdOfCovMat.singularValues();
    Eigen::Matrix2d v = svdOfCovMat.matrixV();
    
    rot = v * u.transpose();

    if (rot.determinant() < 0)
    {
      if (s(0,0) < s(1,0))
        u.col(0) = -u.col(0);
      else
        u.col(1) = -u.col(1);
      rot = v * u.transpose();
    }

    const int idx = 4*t;
    R(idx)   = rot(0,0);
    R(idx+1) = rot(1,0);
    R(idx+2) = rot(0,1);
    R(idx+3) = rot(1,1);
  }
}

void LGARAP_LIMSolver2D::debugOutput(std::stringstream& info)
{
  std::cout << "AR:" << info.str() << "\n"; 
}

void LGARAP_LIMSolver2D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numVertices = mesh->InitalVertices->rows();
  const int numTriangles = mesh->Triangles->rows();

  CotanWeights.resize(numTriangles,3);
  igl::cotmatrix_entries(*mesh->InitalVertices,*mesh->Triangles,CotanWeights);

  // Create matrices L, K
  Eigen::SparseMatrix<double> B, tempL, tempK, restV;
  tempL.resize(numVertices,numVertices);
  tempK.resize(2*numTriangles,numVertices);
  std::vector<Eigen::Triplet<double> > BTriplets, KTriplets, VTriplets, LTriplets;
  for(int t=0;t<numTriangles;t++)
  {
    B.resize(numVertices,3);
    restV.resize(numVertices,2);
    BTriplets.clear();
    VTriplets.clear();
    for(int i=0;i<3;i++)
    {
      const int idx0 = TriEdgeVertices[i][0];
      const int idx1 = TriEdgeVertices[i][1];
      const int vIdx0 = mesh->Triangles->coeff(t,idx0);
      const int vIdx1 = mesh->Triangles->coeff(t,idx1);

      // Create incident matrix B_i for i'th triangle
      BTriplets.push_back(Eigen::Triplet<double>(vIdx0,i,1));
      BTriplets.push_back(Eigen::Triplet<double>(vIdx1,i,-1));

      // Create 3D tet rest pose vertex matrix
      Eigen::Vector2d v0, v1;
      if(mesh->IsCorotatedTriangles)
      {
        v0 = mesh->CorotatedTriangles->block<1,2>(t,2*idx0);
        v1 = mesh->CorotatedTriangles->block<1,2>(t,2*idx1);
      }
      else
      {
        v0 = mesh->InitalVertices->row(mesh->Triangles->coeff(t,idx0)).block<1,2>(0,0);
        v1 = mesh->InitalVertices->row(mesh->Triangles->coeff(t,idx1)).block<1,2>(0,0);
      }
      VTriplets.push_back(Eigen::Triplet<double>(vIdx0,0,v0(0)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx0,1,v0(1)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx1,0,v1(0)));
      VTriplets.push_back(Eigen::Triplet<double>(vIdx1,1,v1(1)));
    }
    B.setFromTriplets(BTriplets.begin(),BTriplets.end());
    restV.setFromTriplets(VTriplets.begin(),VTriplets.end());

    Eigen::SparseMatrix<double> Cm(3,3);
    Cm.insert(0,0) = CotanWeights(t,0);
    Cm.insert(1,1) = CotanWeights(t,1);
    Cm.insert(2,2) = CotanWeights(t,2);

    // Create B*C*B'
    Eigen::SparseMatrix<double> BCBT = B*Cm*B.transpose();
    
    // Stack up K temp matrix
    Eigen::SparseMatrix<double> tempK = 0.5*restV.transpose()*BCBT;
    for (int k=0;k<tempK.outerSize();++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(tempK,k);it;++it)
      {
        KTriplets.push_back(Eigen::Triplet<double>(2*t+it.row(),it.col(),it.value()));
      }
    }

    // Sum up L temp matrix
    for (int k=0;k<BCBT.outerSize();++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(BCBT,k);it;++it)
      {
        int row = 2*it.row();
        int col = 2*it.col();
        LTriplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
        LTriplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
      }
    }
  }

  // Create L matrix
  L.resize(numVariables,numVariables);
  L.setFromTriplets(LTriplets.begin(),LTriplets.end());

  // Create K matrix
  tempK.setFromTriplets(KTriplets.begin(),KTriplets.end());

  K.resize(4*numTriangles,2*numVertices);
  KTriplets.clear();
  for (int k=0;k<tempK.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(tempK,k);it;++it)
    {
      const int row = 2*it.row();
      const int col = 2*it.col();
      KTriplets.push_back(Eigen::Triplet<double>(row,col,it.value()));
      KTriplets.push_back(Eigen::Triplet<double>(row+1,col+1,it.value()));
    }
  }
  K.setFromTriplets(KTriplets.begin(), KTriplets.end());

  RestPoseEdges.resize(3*numTriangles,3); 
  Eigen::RowVector2d p[3];
  for (int t=0;t<numTriangles;t++)
  {
    if(mesh->IsCorotatedTriangles)
    {
      p[0] = mesh->CorotatedTriangles->block<1,2>(t,0);
      p[1] = mesh->CorotatedTriangles->block<1,2>(t,2);
      p[2] = mesh->CorotatedTriangles->block<1,2>(t,4);
    }
    else
    {
      p[0] = mesh->InitalVertices->row(mesh->Triangles->coeff(t,0)).block<1,2>(0,0);
      p[1] = mesh->InitalVertices->row(mesh->Triangles->coeff(t,1)).block<1,2>(0,0);
      p[2] = mesh->InitalVertices->row(mesh->Triangles->coeff(t,2)).block<1,2>(0,0);
    }

    for(int i=0;i<3;i++)
      RestPoseEdges.block<2,1>(2*t,i) = p[TriEdgeVertices[i][1]]-p[TriEdgeVertices[i][0]];
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
  int numVariables = numVertices*2;
  Eigen::Matrix<double,Eigen::Dynamic,1> restPose(numVariables);
  for(int n=0;n<numVertices;n++)
  {
    for(int i=0;i<2;i++)
      restPose[n*2+i] = mesh->InitalVertices->coeff(n,i);
  }
  constantEnergyPart = 0.5*restPose.transpose()*L*restPose;
}

double LGARAP_LIMSolver2D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // ARAP global step energy
  double xTLx = x.transpose()*L*x;
  double KTRx = R.transpose()*K*x;
  return 0.5*xTLx - KTRx + constantEnergyPart;
}

void LGARAP_LIMSolver2D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // ARAP global step
  grad = L*x - K.transpose()*R;
}

void LGARAP_LIMSolver2D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
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
