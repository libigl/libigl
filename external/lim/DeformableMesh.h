// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#pragma once

#ifndef DEFORMABLE_MESH_H
#define DEFORMABLE_MESH_H

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class DeformableMesh
{
public:

  // The following data is needed to run the solver

  // vertex positions of initial rest pose
  Eigen::Matrix<double,Eigen::Dynamic,3>* InitalVertices;
  
  // vertex positions of current mesh
  Eigen::Matrix<double,Eigen::Dynamic,3>* DeformedVertices;
  
  // new vertex positions after a solve
  Eigen::Matrix<double,Eigen::Dynamic,3>* PredictedVertices;
  
  // indices of all border vertices
  Eigen::Matrix<int,Eigen::Dynamic,1>* BorderVertices;
  
  // positional constraint data: Cv-d
  // C: linear positional constraint matrix of dimension: (#constraints)x(v*dim)
  Eigen::SparseMatrix<double>* ConstraintMatrix;
  // d: target positions vector of dimension #constraints
  Eigen::Matrix<double,Eigen::Dynamic,1>* ConstraintTargets;
  
  double Radius;
  Eigen::Vector3f Center;
  Eigen::Vector3f BoxSize;
  Eigen::Vector3f BoxCornerMin;
  Eigen::Vector3f BoxCornerMax;

  double NodeSize;
  double EPS1;
  double EPS3;

  DeformableMesh() {};
  virtual ~DeformableMesh() {};

  virtual void InitMesh() = 0;
  virtual void UpdateMesh() = 0;
};

#endif