// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "GreenStrain_LIMSolver3D.h"
#include "TetrahedronMesh.h"

GreenStrain_LIMSolver3D::GreenStrain_LIMSolver3D()
{
  Beta = 100;
}

GreenStrain_LIMSolver3D::~GreenStrain_LIMSolver3D()
{
}

void GreenStrain_LIMSolver3D::debugOutput(std::stringstream& info)
{
  std::cout << "GS:" << info.str() << "\n"; 
}

void GreenStrain_LIMSolver3D::prepareProblemData(std::vector<int>& hessRowIdx, std::vector<int>& hessColIdx)
{
  const int numNodes = mesh->InitalVertices->rows();

  // Compute deformation gradients
  int numTets = mesh->Tetrahedra->rows();
  Ms.resize(4,3*numTets);
  MMTs.resize(4,4*numTets);

  Eigen::Matrix<double,4,3> SelectorM;
  SelectorM.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
  SelectorM.row(3) = Eigen::Vector3d::Ones()*-1;
  
  for(int t=0;t<numTets;t++)
  {
    Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);

    Eigen::Vector3d A = mesh->InitalVertices->row(mesh->Tetrahedra->coeff(t,0)).cast<double>();
    Eigen::Vector3d B = mesh->InitalVertices->row(mesh->Tetrahedra->coeff(t,1)).cast<double>();
    Eigen::Vector3d C = mesh->InitalVertices->row(mesh->Tetrahedra->coeff(t,2)).cast<double>();
    Eigen::Vector3d D = mesh->InitalVertices->row(mesh->Tetrahedra->coeff(t,3)).cast<double>();

    Eigen::Matrix3d V;
    V << A-D,B-D,C-D;
    
    Eigen::Matrix<double,4,3> Mtemp = SelectorM*V.inverse().cast<double>();
    Ms.block<4,3>(0,3*t) = Mtemp;
    MMTs.block<4,4>(0,4*t) = Mtemp*Mtemp.transpose();
  }
}

double GreenStrain_LIMSolver3D::computeFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  // green strain energy
  double shape = 0;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  for(int t=0;t<mesh->Tetrahedra->rows();t++)
  {
    Eigen::Vector3d A(x[TetrahedronVertexIdx.coeff(0,t)],x[TetrahedronVertexIdx.coeff(1,t)],x[TetrahedronVertexIdx.coeff(2,t)]);
    Eigen::Vector3d B(x[TetrahedronVertexIdx.coeff(3,t)],x[TetrahedronVertexIdx.coeff(4,t)],x[TetrahedronVertexIdx.coeff(5,t)]);
    Eigen::Vector3d C(x[TetrahedronVertexIdx.coeff(6,t)],x[TetrahedronVertexIdx.coeff(7,t)],x[TetrahedronVertexIdx.coeff(8,t)]);
    Eigen::Vector3d D(x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(11,t)]);

    Eigen::Matrix<double,3,4> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;
    V.col(3) = D;

    Eigen::Matrix3d F = V*Ms.block<4,3>(0,3*t);
    Eigen::Matrix3d E = (F.transpose()*F - I);
    shape += E.squaredNorm()*Divider;
  }

  return shape;
}

void GreenStrain_LIMSolver3D::computeGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  // green strain energy
  for(int t=0;t<mesh->Tetrahedra->rows();t++)
  {
    Eigen::Vector3d A(x[TetrahedronVertexIdx.coeff(0,t)],x[TetrahedronVertexIdx.coeff(1,t)],x[TetrahedronVertexIdx.coeff(2,t)]);
    Eigen::Vector3d B(x[TetrahedronVertexIdx.coeff(3,t)],x[TetrahedronVertexIdx.coeff(4,t)],x[TetrahedronVertexIdx.coeff(5,t)]);
    Eigen::Vector3d C(x[TetrahedronVertexIdx.coeff(6,t)],x[TetrahedronVertexIdx.coeff(7,t)],x[TetrahedronVertexIdx.coeff(8,t)]);
    Eigen::Vector3d D(x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(11,t)]);

    Eigen::Matrix<double,3,4> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;
    V.col(3) = D;

    // jacobian(E) = 4(VMM'V'VMM' - VMM')
    Eigen::Matrix<double,3,4> VMMT = V*MMTs.block<4,4>(0,4*t);
    Eigen::Matrix<double,3,4> T = 4*(VMMT*V.transpose()*VMMT - VMMT);

    for(int i=0;i<12;i++)
      grad[TetrahedronVertexIdx.coeff(i,t)] += T.coeff(i)*Divider;
  }
}

void GreenStrain_LIMSolver3D::computeHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double*,Eigen::Dynamic,1>& hess)
{
  // green strain tensor energy
  Eigen::Matrix<double,3,4> S;
  for(int t=0;t<mesh->Tetrahedra->rows();t++)
  {
    Eigen::Vector3d A(x[TetrahedronVertexIdx.coeff(0,t)],x[TetrahedronVertexIdx.coeff(1,t)],x[TetrahedronVertexIdx.coeff(2,t)]);
    Eigen::Vector3d B(x[TetrahedronVertexIdx.coeff(3,t)],x[TetrahedronVertexIdx.coeff(4,t)],x[TetrahedronVertexIdx.coeff(5,t)]);
    Eigen::Vector3d C(x[TetrahedronVertexIdx.coeff(6,t)],x[TetrahedronVertexIdx.coeff(7,t)],x[TetrahedronVertexIdx.coeff(8,t)]);
    Eigen::Vector3d D(x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(11,t)]);

    Eigen::Matrix<double,3,4> V;
    V.col(0) = A;
    V.col(1) = B;
    V.col(2) = C;
    V.col(3) = D;

    // hessian(E) = 4*r_x'*((SMM'V'V+VMM'*(V'S+SV))*MM' - SMM')*c_x
    Eigen::Matrix<double,4,4> VTV = V.transpose()*V;
    Eigen::Matrix<double,4,4> MMT = MMTs.block<4,4>(0,4*t);
    Eigen::Matrix<double,3,4> VMMT = V*MMT;
    Eigen::Matrix<double,4,4> MMTVTV = MMT*VTV;

    int numElem = 0;
    for(int r=0;r<12;r++)
    {
      S = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(3,4);
      S.coeffRef(r) = 1;

      Eigen::Matrix<double,3,4> Temp = 4*((S*MMTVTV + VMMT*(V.transpose()*S+S.transpose()*V))*MMT - S*MMT);
      
      for(int c=r;c<12;c++)
        *denseHessianCoeffs(numElem++,t) += Temp.coeff(c)*Divider;
    }
  }
}
