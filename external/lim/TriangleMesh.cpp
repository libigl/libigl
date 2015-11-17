// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "TriangleMesh.h"

TriangleMesh::TriangleMesh()
{
  InitalVertices = NULL;
  DeformedVertices = NULL;
  BorderVertices = NULL;
  ConstraintMatrix = NULL;
  ConstraintTargets = NULL;

  IsCorotatedTriangles = false;

  Radius = 0;
  Center = Eigen::Vector3f(0,0,0);
}

TriangleMesh::~TriangleMesh()
{
}
  
void TriangleMesh::InitMesh()
{
  UpdateMesh();

  // Find min triangle area
  double minArea = std::numeric_limits<double>::max();
  double avgArea = 0;
  for(int t=0;t<Triangles->rows();t++)
  {
    Eigen::Vector3d A = InitalVertices->row(Triangles->coeff(t,0)).cast<double>();
    Eigen::Vector3d B = InitalVertices->row(Triangles->coeff(t,1)).cast<double>();
    Eigen::Vector3d C = InitalVertices->row(Triangles->coeff(t,2)).cast<double>();

    double area = ((A-C).cross(B-C)).norm();

    if(area <= 0)
      std::cerr << "element " << t << " has zero or negative area!\n";

    avgArea += area;

    if(area < minArea)
      minArea = area;
  }
  avgArea /= Triangles->rows();

  EPS1 = 10e-5;
  EPS3 = minArea*EPS1;
}

void TriangleMesh::UpdateMesh()
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
