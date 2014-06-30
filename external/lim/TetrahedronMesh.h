// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef TETRAHEDRON_OBJECT_H
#define TETRAHEDRON_OBJECT_H

#include "DeformableMesh.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class TetrahedronMesh : public DeformableMesh
{
public:
  Eigen::Matrix<int,Eigen::Dynamic,4>* Tetrahedra;
  
  TetrahedronMesh();
  ~TetrahedronMesh();
  
  void InitMesh();
  void UpdateMesh();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
};

#endif