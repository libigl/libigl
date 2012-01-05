#ifndef IGL_MIN_QUAD_DENSE_H
#define IGL_MIN_QUAD_DENSE_H

#include <Eigen/Core>
#include <Eigen/Dense>

//// debug
//#include <matlabinterface.h>
//Engine *g_pEngine;

#include <EPS.h>

namespace igl
{
	// MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C
	// subject to linear constraints Aeq*Z = Beq
	//
	// Templates:
	//   T  should be a eigen matrix primitive type like float or double
	// Inputs:
	//   A  n by n matrix of quadratic coefficients
	//   B  n by 1 column of linear coefficients
	//   Aeq  m by n list of linear equality constraint coefficients
	//   Beq  m by 1 list of linear equality constraint constant values
	// Outputs:
	//   S  n by (n + m) "solve" matrix, such that S*[B', Beq'] is a solution
	// Returns true on success, false on error
	template <typename T>
	void min_quad_dense_precompute(
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Aeq,		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& S)
	{
		typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
                // This threshold seems to matter a lot but I'm not sure how to
                // set it
		const T treshold = igl::FLOAT_EPS;
		//const T treshold = igl::DOUBLE_EPS;

		const int n = A.rows();
		assert(A.cols() == n);
		const int m = Aeq.rows();
		assert(Aeq.cols() == n);

		// Lagrange multipliers method:
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> LM(n + m, n + m);
		LM.block(0, 0, n, n) = A;
		LM.block(0, n, n, m) = Aeq.transpose();
		LM.block(n, 0, m, n) = Aeq;
		LM.block(n, n, m, m).setZero();

		typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec; 
		Vec singValues;
		Eigen::JacobiSVD<Mat> svd;
		svd.compute(LM, Eigen::ComputeFullU | Eigen::ComputeFullV );
		const Mat& u = svd.matrixU();
		const Mat& v = svd.matrixV();
		const Vec& singVals = svd.singularValues();

		Vec pi_singVals(n + m);
		int zeroed = 0;
		for (int i=0; i<n + m; i++)
		{
			T sv = singVals(i, 0);
			assert(sv >= 0);			
                        //printf("sv: %lg ? %lg\n",(double) sv,(double)treshold);
			if (sv > treshold) pi_singVals(i, 0) = T(1) / sv;
			else 
			{
				pi_singVals(i, 0) = T(0);
				zeroed++;
			}
		}

		printf("min_quad_dense_precompute: %i singular values zeroed\n", zeroed);
		Eigen::DiagonalMatrix<T, Eigen::Dynamic> pi_diag(pi_singVals);

		Mat LMpinv = v * pi_diag * u.transpose();
		S = LMpinv.block(0, 0, n, n + m);

		//// debug:
		//mlinit(&g_pEngine);
		//
		//mlsetmatrix(&g_pEngine, "A", A);
		//mlsetmatrix(&g_pEngine, "Aeq", Aeq);
		//mlsetmatrix(&g_pEngine, "LM", LM);
		//mlsetmatrix(&g_pEngine, "u", u);
		//mlsetmatrix(&g_pEngine, "v", v);
		//MatrixXd svMat = singVals;
		//mlsetmatrix(&g_pEngine, "singVals", svMat);
		//mlsetmatrix(&g_pEngine, "LMpinv", LMpinv);
		//mlsetmatrix(&g_pEngine, "S", S);

		//int hu = 1;
	}
}

#endif
