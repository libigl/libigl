// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "GreenStrain_LIMSolver2D.h"
#include "TriangleMesh.h"

GreenStrain_LIMSolver2D::GreenStrain_LIMSolver2D()
{
  Beta = 100;
}

GreenStrain_LIMSolver2D::~GreenStrain_LIMSolver2D()
{
}

void GreenStrain_LIMSolver2D::debugOutput(std::stringstream& info)
{
  std::cout << "GS:" << info.str() << "\n"; 
}

void GreenStrain_LIMSolver2D::getProblemSize()
{
  numVariables = mesh->InitalVertices->rows()*2;
}

void GreenStrain_LIMSolver2D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  // Compute deformation gradients
  int numTets = mesh->Triangles->rows();
  Ms.resize(3,2*numTets);
  MMTs.resize(3,3*numTets);

  Eigen::Matrix<double,3,2> SelectorM;
  SelectorM.block<2,2>(0,0) = Eigen::Matrix<double,2,2>::Identity();
  SelectorM.row(2) = Eigen::Vector2d::Ones()*-1;
  
  for(int t=0;t<numTets;t++)
  {
    Eigen::Vector2d A,B,C;
    if(mesh->IsCorotatedTriangles)
    {
      A = mesh->CorotatedTriangles->row(t).block<1,2>(0,0).cast<double>();
      B = mesh->CorotatedTriangles->row(t).block<1,2>(0,2).cast<double>();
      C = mesh->CorotatedTriangles->row(t).block<1,2>(0,4).cast<double>();	
    }
    else
    {
      A = mesh->InitalVertices->row(mesh->Triangles->coeff(t,0)).block<1,2>(0,0).cast<double>();
      B = mesh->InitalVertices->row(mesh->Triangles->coeff(t,1)).block<1,2>(0,0).cast<double>();
      C = mesh->InitalVertices->row(mesh->Triangles->coeff(t,2)).block<1,2>(0,0).cast<double>();
    }

    Eigen::Matrix2d V;
    V << A-C,B-C;
    
    Eigen::Matrix<double,3,2> Mtemp = SelectorM*V.inverse().cast<double>();
    Ms.block<3,2>(0,2*t) = Mtemp;
    MMTs.block<3,3>(0,3*t) = Mtemp*Mtemp.transpose();
  }
}

double GreenStrain_LIMSolver2D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // green strain energy
  double shape = 0;
  Eigen::Matrix<double,2,2> I = Eigen::Matrix<double,2,2>::Identity();
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    Eigen::Vector2d A(x[TriangleVertexIdx.coeff(0,t)],x[TriangleVertexIdx.coeff(1,t)]);
    Eigen::Vector2d B(x[TriangleVertexIdx.coeff(2,t)],x[TriangleVertexIdx.coeff(3,t)]);
    Eigen::Vector2d C(x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(5,t)]);

    Eigen::Matrix<double,2,3> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;

    Eigen::Matrix<double,2,2> F = V*Ms.block<3,2>(0,2*t);
    Eigen::Matrix<double,2,2> E = (F.transpose()*F - I);
    shape += E.squaredNorm()*Divider;
  }

    return shape;
}

void GreenStrain_LIMSolver2D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // green strain energy
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    Eigen::Vector2d A(x[TriangleVertexIdx.coeff(0,t)],x[TriangleVertexIdx.coeff(1,t)]);
    Eigen::Vector2d B(x[TriangleVertexIdx.coeff(2,t)],x[TriangleVertexIdx.coeff(3,t)]);
    Eigen::Vector2d C(x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(5,t)]);

    Eigen::Matrix<double,2,3> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;

    // jacobian(E) = 4(VMM'V'VMM' - VMM')
    Eigen::Matrix<double,2,3> VMMT = V*MMTs.block<3,3>(0,3*t);
    Eigen::Matrix<double,2,3> T = 4*(VMMT*V.transpose()*VMMT - VMMT);

    for(int i=0;i<6;i++)
      grad[TriangleVertexIdx.coeff(i,t)] += T.coeff(i)*Divider;
  }
}

void GreenStrain_LIMSolver2D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // green strain tensor energy
  Eigen::Matrix<double,2,3> S;
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    Eigen::Vector2d A(x[TriangleVertexIdx.coeff(0,t)],x[TriangleVertexIdx.coeff(1,t)]);
    Eigen::Vector2d B(x[TriangleVertexIdx.coeff(2,t)],x[TriangleVertexIdx.coeff(3,t)]);
    Eigen::Vector2d C(x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(5,t)]);

    Eigen::Matrix<double,2,3> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;

    // hessian(E) = 4*r_x'*((SMM'V'V+VMM'*(V'S+SV))*MM' - SMM')*c_x
    Eigen::Matrix3d VTV = V.transpose()*V;
    Eigen::Matrix3d MMT = MMTs.block<3,3>(0,3*t);
    Eigen::Matrix<double,2,3> VMMT = V*MMT;
    Eigen::Matrix3d MMTVTV = MMT*VTV;

    int numElem = 0;
    for(int r=0;r<6;r++)
    {
      S = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(2,3);
      S.coeffRef(r) = 1;

      Eigen::Matrix<double,2,3> Temp = 4*((S*MMTVTV + VMMT*(V.transpose()*S+S.transpose()*V))*MMT - S*MMT);
      
      for(int c=r;c<6;c++)
        *denseHessianCoeffs(numElem++,t) += Temp.coeff(c)*Divider;
    }
  }
}
