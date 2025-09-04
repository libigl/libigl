// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "shortest_edge_and_midpoint.h"
#include <iostream>
#include <Eigen/LU>

IGL_INLINE void igl::shortest_edge_and_midpoint(
  const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)
{
  cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
  p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
}



IGL_INLINE void igl::shortest_edge_and_midpoint3(
	const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F/*F*/,
	const Eigen::MatrixXi& E,
	const Eigen::VectorXi& EMAP/*EMAP*/,
	const Eigen::MatrixXi& EF/*EF*/,
	const Eigen::MatrixXi& EI/*EI*/,
	const Eigen::MatrixXd& PD2,
	const Eigen::VectorXd& PV2,
	const std::vector<Eigen::Matrix4d>& q_matrices,
	const double& avg_edges_lenth,
	double& cost,
	Eigen::RowVectorXd& p)
{
	const int v1_idx = E(e, 0);
	const int v2_idx = E(e, 1);
	double base_cost = (V.row(E(e, 0)) - V.row(E(e, 1))).norm();
	// --- 1. QEM ���� ---
	// �ϲ�Q����
	Eigen::Matrix4d Q_new = q_matrices[v1_idx] + q_matrices[v2_idx];
	Eigen::Matrix4d Q_pre = Q_new;
	Eigen::Vector4d b(0.0, 0.0, 0.0, 1.0);
	Q_pre(3, 0) = 0.0;
	Q_pre(3, 1) = 0.0;
	Q_pre(3, 2) = 0.0;
	Q_pre(3, 3) = 1.0;
	Eigen::FullPivLU<Eigen::Matrix4d> lu(Q_pre);
	Eigen::Vector4d p_homo;
	// if is invertible, solve the linear equation
	if (lu.isInvertible())
	{
		p_homo = Q_pre.inverse() * b;
	}
	// else take the midpoint
	else
	{
		p = 0.5 * (V.row(v1_idx) + V.row(v2_idx));
		p_homo << p.transpose(), 1.0;
	}
	//// ʹ�� Eigen �����������
	//Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	//// �������Ƿ���� (ͨ���������ֵ)
	//if (svd.singularValues().minCoeff() > 1e-6) {
	//	p = svd.solve(b).transpose();
	//}
	//else {
	//	// ��������棬�˻ص��е�
	//	p = 0.5 * (V.row(v1_idx) + V.row(v2_idx));
	//}
	assert(fabs(p_homo(3) - 1) < 1e-6);
	// ����QEM����
	p = p_homo.head<3>().transpose();

	double cost_qem = p_homo.transpose() * Q_new * p_homo;

	// --- 2. ���ʺͱ߽�ͷ� (��֮ǰ�汾��ͬ) ---
	// ... (ʡ�ԣ���������һ�汾��ȫ��ͬ) ...o
	const Eigen::RowVector3d edge_dir = (V.row(v2_idx) - V.row(v1_idx)).normalized();
	Eigen::RowVector3d d_min1 = PD2.row(v1_idx);
	Eigen::RowVector3d d_min2 = PD2.row(v2_idx);
	const double v_min1 = fabs(PV2(v1_idx));
	const double v_min2 = fabs(PV2(v2_idx));
	double curvature_penalty = fabs(std::max(1 - std::abs(edge_dir.dot(d_min1)), 1 - std::abs(edge_dir.dot(d_min2))));
	/*if (v_min1 < 1e-6 && v_min2<1e-6) curvature_penalty = 0;*/
	
	
	double boundary_penalty = (EF(e, 0) == -1 || EF(e, 1) == -1) ? 1000.0 : 0.0;
	const double w_curve = 1;

	// --- 3. ������մ��� ---
	/*cost = cost_qem * (w_curve * curvature_penalty) + boundary_penalty;*/
	/*cost = (cost_qem + 0.1 * base_cost * base_cost) * (1+w_curve * curvature_penalty) + boundary_penalty;*/
	cost = cost_qem + w_curve* curvature_penalty;
}