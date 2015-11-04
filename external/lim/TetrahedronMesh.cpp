// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "TetrahedronMesh.h"

TetrahedronMesh::TetrahedronMesh()
{
  InitalVertices = NULL;
  DeformedVertices = NULL;
  BorderVertices = NULL;
  ConstraintMatrix = NULL;
  ConstraintTargets= NULL;

  Radius = 0;
  Center = Eigen::Vector3f(0,0,0);
}

TetrahedronMesh::~TetrahedronMesh()
{
}
  
void TetrahedronMesh::InitMesh()
{
  UpdateMesh();

  // Find min tet volume
  double minVol = std::numeric_limits<double>::max();
  for(int t=0;t<Tetrahedra->rows();t++)
  {
    Eigen::Vector3d A = InitalVertices->row(Tetrahedra->coeff(t,0)).cast<double>();
    Eigen::Vector3d B = InitalVertices->row(Tetrahedra->coeff(t,1)).cast<double>();
    Eigen::Vector3d C = InitalVertices->row(Tetrahedra->coeff(t,2)).cast<double>();
    Eigen::Vector3d D = InitalVertices->row(Tetrahedra->coeff(t,3)).cast<double>();

    Eigen::Vector3d a = A-D;
    Eigen::Vector3d b = B-D;
    Eigen::Vector3d c = C-D;

    double vol = a.dot(c.cross(b));

    if(vol <= 0)
      std::cerr << "element " << t << " has zero or negative volume!\n";

    if(vol < minVol)
      minVol = vol;
  }

  EPS1 = 10e-5;
  EPS3 = minVol*EPS1;
}

void TetrahedronMesh::UpdateMesh()
{
  Eigen::Vector3f minVal = DeformedVertices->row(0).cast<float>();
  Eigen::Vector3f maxVal = DeformedVertices->row(0).cast<float>();
  for(int i=0;i<DeformedVertices->rows();i++)
  {
    Eigen::Vector3f v = DeformedVertices->row(i).cast<float>();

    for(int i=0;i<3;i++)
    {
      if(v(i) < minVal(i))
        minVal(i) = v(i);

      if(v(i) > maxVal(i))
        maxVal(i) = v(i);	
    }
  }

  BoxCornerMin = minVal;
  BoxCornerMax = maxVal;
  BoxSize = maxVal -minVal;
  Radius = BoxSize.norm()*0.5f;
  Center = ((maxVal+minVal)*0.5f).cast<float>();
  NodeSize = Radius*0.02;
}
