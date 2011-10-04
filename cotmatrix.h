//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef COTMATRIX_H
#define COTMATRIX_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace igl 
{
    void computeCotWeights(Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double& cotAlpha, double& cotBeta, double& cotGamma)
    {
        Eigen::Vector3d v12 = v2-v1;
        Eigen::Vector3d v13 = v3-v1;
        Eigen::Vector3d v23 = v3-v2;
        
        double halfArea = (v12.cross(v13)).norm();//squaredNorm();
        
        //improve numerical stability
        const double cotTolerance = 1e-10;
        if(halfArea < cotTolerance)
        {
            std::cout << "Cot weights are close to singular!" << std::endl;
            halfArea = cotTolerance;
        }
        
        cotAlpha = (v12.dot(v13)) / halfArea /2;
        cotBeta =  (v23.dot(-v12)) / halfArea /2;
        cotGamma = (-v23.dot(-v13)) / halfArea /2;
    }
    
    void cotmatrix (Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::SparseMatrix<double>& L, bool cotan = true)
    {
        Eigen::DynamicSparseMatrix<double, Eigen::RowMajor> dyn_L (V.rows(), V.rows());
        
        for (unsigned i = 0; i < F.rows(); i++)
        {
            int vi1 = F(i,0);
            int vi2 = F(i,1);
            int vi3 = F(i,2);
            
            if (cotan)
            {
                Eigen::Vector3d v1 (V(vi1,0), V(vi1,1), V(vi1,2));
                Eigen::Vector3d v2 (V(vi2,0), V(vi2,1), V(vi2,2));
                Eigen::Vector3d v3 (V(vi3,0), V(vi3,1), V(vi3,2));
                double cot_a, cot_b, cot_g;
                computeCotWeights (v1, v2, v3, cot_a, cot_b, cot_g);
                
                dyn_L.coeffRef (vi1, vi2) += cot_g;
                dyn_L.coeffRef (vi2, vi1) += cot_g;
                dyn_L.coeffRef (vi2, vi3) += cot_a;
                dyn_L.coeffRef (vi3, vi2) += cot_a;
                dyn_L.coeffRef (vi3, vi1) += cot_b;
                dyn_L.coeffRef (vi1, vi3) += cot_b;
            }
            else
            {
                dyn_L.coeffRef (vi1, vi2) += 1.0;
                dyn_L.coeffRef (vi2, vi3) += 1.0;
                dyn_L.coeffRef (vi3, vi1) += 1.0;
            }
        }
        for (int k=0; k < dyn_L.outerSize(); ++k)
        {
            double tmp = 0.0f;
            for (Eigen::DynamicSparseMatrix<double, Eigen::RowMajor>::InnerIterator it (dyn_L, k); it; ++it)
                tmp += it.value ();
            dyn_L.coeffRef (k,k) = -tmp;
        }
        L = Eigen::SparseMatrix<double> (dyn_L);
    }
}

#endif