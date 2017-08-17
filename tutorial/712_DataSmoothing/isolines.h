
#include <vector>
#include <limits>
#include <stdlib.h>

#include <igl/remove_unreferenced.h>


static void isolines(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXd& z,
    const int grads,
    Eigen::MatrixXd& isoV,
    Eigen::MatrixXi& isoE)
{
    const double min = z.minCoeff(), max = z.maxCoeff();

    //Following http://www.alecjacobson.com/weblog/?p=2529
    Eigen::VectorXd iso(grads+1);
    for(int i=0; i<iso.size(); ++i)
        iso(i) = double(i)/double(grads)*(max-min) + min;

    Eigen::MatrixXd t12(F.rows(), iso.size()), t23(F.rows(), iso.size()),
        t31(F.rows(), iso.size());
    for(int i=0; i<t12.rows(); ++i) {
        const double z1=z(F(i,0)), z2=z(F(i,1)), z3=z(F(i,2));
        const double s12 = z2-z1;
        const double s23 = z3-z2;
        const double s31 = z1-z3;
        for(int j=0; j<t12.cols(); ++j) {
            t12(i,j) = (iso(j)-z1) / s12;
            t23(i,j) = (iso(j)-z2) / s23;
            t31(i,j) = (iso(j)-z3) / s31;
            if(t12(i,j)<0 || t12(i,j)>1)
                t12(i,j) = std::numeric_limits<double>::quiet_NaN();
            if(t23(i,j)<0 || t23(i,j)>1)
                t23(i,j) = std::numeric_limits<double>::quiet_NaN();
            if(t31(i,j)<0 || t31(i,j)>1)
                t31(i,j) = std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::vector<int> F12, F23, F31, I12, I23, I31;
    for(int i=0; i<t12.rows(); ++i) {
        for(int j=0; j<t12.cols(); ++j) {
            if(std::isfinite(t23(i,j)) && std::isfinite(t31(i,j))) {
                F12.push_back(i);
                I12.push_back(j);
            }
            if(std::isfinite(t31(i,j)) && std::isfinite(t12(i,j))) {
                F23.push_back(i);
                I23.push_back(j);
            }
            if(std::isfinite(t12(i,j)) && std::isfinite(t23(i,j))) {
                F31.push_back(i);
                I31.push_back(j);
            }
        }
    }

    const int K = F12.size()+F23.size()+F31.size();
    isoV.resize(2*K, 3);
    int b = 0;
    for(int i=0; i<F12.size(); ++i) {
        isoV.row(b+i) = (1.-t23(F12[i],I12[i]))*V.row(F(F12[i],1))
            + t23(F12[i],I12[i])*V.row(F(F12[i],2));
        isoV.row(K+b+i) = (1.-t31(F12[i],I12[i]))*V.row(F(F12[i],2))
            + t31(F12[i],I12[i])*V.row(F(F12[i],0));
    }
    b += F12.size();
    for(int i=0; i<F23.size(); ++i) {
        isoV.row(b+i) = (1.-t31(F23[i],I23[i]))*V.row(F(F23[i],2))
            + t31(F23[i],I23[i])*V.row(F(F23[i],0));
        isoV.row(K+b+i) = (1.-t12(F23[i],I23[i]))*V.row(F(F23[i],0))
            + t12(F23[i],I23[i])*V.row(F(F23[i],1));
    }
    b += F23.size();
    for(int i=0; i<F31.size(); ++i) {
        isoV.row(b+i) = (1.-t12(F31[i],I31[i]))*V.row(F(F31[i],0))
            + t12(F31[i],I31[i])*V.row(F(F31[i],1));
        isoV.row(K+b+i) = (1.-t23(F31[i],I31[i]))*V.row(F(F31[i],1))
            + t23(F31[i],I31[i])*V.row(F(F31[i],2));
    }

    isoE.resize(K,2);
    for(int i=0; i<K; ++i)
        isoE.row(i) << i, K+i;

}

