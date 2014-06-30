// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef TRIANGLE_OBJECT_H
#define TRIANGLE_OBJECT_H

#include "DeformableMesh.h"

#include <Eigen/Dense>

class TriangleMesh : public DeformableMesh
{
public:
  // vertex indices of each triangle
  Eigen::Matrix<int,Eigen::Dynamic,3>* Triangles;
  
  // used for parametrization only
  bool IsCorotatedTriangles;
  Eigen::Matrix<double,Eigen::Dynamic,6>* CorotatedTriangles;

  TriangleMesh();
  ~TriangleMesh();
  
  void InitMesh();
  void UpdateMesh();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
};

#endif