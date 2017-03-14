
#include "discrete_exterior_calculus.h"
#include "halfedge_flip.h"
#include "doublearea.h"
#include "cotmatrix_entries.h"

#include <vector>
#include <iostream>
namespace igl
{
    namespace dec2d
    {
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void d0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D)
        {
            using namespace std;
            if (F.cols() != 3)
            {
                cerr << "discrete_exterior_calculus.h: Error: Simplex size ("<<
                    F.cols() << ") not supported" << endl;
                assert(false);
            }

            int m = F.rows();
            int n = V.rows();
            /* build derivative operator matrix
             * for halfedge he = (u, v):
             * D(he, u) = -1, D(he, v) = 1 
             */ 
            typedef Eigen::Triplet<Scalar> T;
            vector<T> d0_coeff;
            d0_coeff.clear();
            d0_coeff.reserve(6 * m);
            for(int i = 0; i < m; i++)
            {
                d0_coeff.push_back(T(3 * i    , F(i, 1), -1));
                d0_coeff.push_back(T(3 * i    , F(i, 2),  1));

                d0_coeff.push_back(T(3 * i + 1, F(i, 2), -1));
                d0_coeff.push_back(T(3 * i + 1, F(i, 0),  1));

                d0_coeff.push_back(T(3 * i + 2, F(i, 0), -1));
                d0_coeff.push_back(T(3 * i + 2, F(i, 1),  1));
            }

            D = Eigen::SparseMatrix<Scalar>(3 * m, n);
            D.setFromTriplets(d0_coeff.begin(), d0_coeff.end());
        }
        
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void d1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D)
        {
            using namespace std;
            if (F.cols() != 3)
            {
                cerr << "discrete_exterior_calculus.h: Error: Simplex size ("<<
                    F.cols() << ") not supported" << endl;
                assert(false);
            }

            int m = F.rows();

            /* get halfedge flip information */
            Eigen::MatrixXi HE;
            halfedge_flip(F, HE);

            /* build direvative operator 
             * if HE is boundary, contribute 1 to facet integration 
             * if not, average the HE on both side and then get facet integration 
             * 
             * t = [u, v, w] 
             * if he0 = (u, v) and he0 has opposite halfedge he1 = (v, u)
             *      D(t, he0) = 0.5, D(t, he1) = -0.5
             * if he0 = (u, v) and he0 is boundary, thus has no opposite halfedge
             *      D(t, he0) = 1.0
             */
            typedef Eigen::Triplet<Scalar> T;
            vector<T> d1_coeff;
            d1_coeff.clear();
            d1_coeff.reserve(6 * m);
            for (int i = 0; i < m; i++)
                for(int j = 0; j < 3; j++)
                {
                    if (HE(i, j) < 0)
                    {
                        d1_coeff.push_back(T(i, 3 * i + j, 1));
                    }else
                    {
                        d1_coeff.push_back(T(i, 3 * i + j, 0.5));
                        d1_coeff.push_back(T(i, HE(i,j),  -0.5));
                    }
                }
            D = Eigen::SparseMatrix<Scalar>(m, 3 * m);
            D.setfromTriplets(d1_coeff.begin(), d1_coeff.end());
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_d0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D)
        {
            d0(V, F, D);
            D = -0.5 * D.transpose();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_d1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D)
        {
            d1(V, F, D);
            D = 2 * D.transpose();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            using namespace std;
            if (F.cols() != 3)
            {
                cerr << "discrete_exterior_calculus.h: Error: Simplex size ("<<
                    F.cols() << ") not supported" << endl;
                assert(false);
            }

            int n = V.rows();
            int m = F.rows();

            Eigen::VectorXd area;
            doublearea(V, F, area);
            area /= 6.0;  /* 2.0 for doubled area, 3.0 for three vertices */
            
            Eigen::VectorXd v_area = Eigen::VectorXd::Zero(n);
            parallel_for(
                    m,
                    [&F, &area, &v_area](const int i)
                    {
                        v_area(F(i, 0)) += area(i);
                        v_area(F(i, 1)) += area(i);
                        v_area(F(i, 2)) += area(i);
                    },
                    1000);
            HS = v_area.asDiagonal();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            using namespace std;
            using namespace Eigen;
            if (F.cols() != 3)
            {
                cerr << "discrete_exterior_calculus.h: Error: Simplex size ("<<
                    F.cols() << ") not supported" << endl;
                assert(false);
            }

            int m = F.rows();
            /* get halfedge flip information */
            Eigen::MatrixXi HE;
            halfedge_flip(F, HE);

            /* get cot value of each halfedge's corresponding angle */
            Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> Cot;
            cotmatrix_entries(V, F, Cot);
            
            /* get transform between halfedges and dual halfedges */
            Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> HS_diag = Cot;
            for(int i = 0; i < m; i++)
                for(int j = 0; j < 3; j++)
                {   /* if halfedge is boundary, contribute cot() to dual halfedge integration 
                     * if not, average the cot() on both side and then get halfedge integration 
                     */
                    if (HE(i, j) < 0) 
                    {
                        HS_diag(i, j) += Cot(i, j);
                    }else
                    {
                        HS_diag(i, j) += Cot(HE(i, j));  /* only when Cot is RowMajor */
                    }
                }
            HS_diag *= 0.5;
            HS_diag.resize(HS_diag.size(), 1);  /* only when HS_diag is RowMajor */
            HS = HS_diag.asDiagonal();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs2(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            using namespace std;
            if (F.cols() != 3)
            {
                cerr << "discrete_exterior_calculus.h: Error: Simplex size ("<<
                    F.cols() << ") not supported" << endl;
                assert(false);
            }
            Eigen::VectorXd area;
            doublearea(V, F, area);
            area /= 2.0;  /* 2.0 for doubled area */
            HS = area.asDiagonal();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            hs0(V, F, HS);
            auto v_area = HS.diagonal();
            v_area = 1.0 / v_area.array();
            HS = v_area.asDiagonal();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            hs1(V, F, HS);
            auto he_length = HS.diagonal();
            he_length = - 1.0 / he_length.array();
            HS = he_length.asDiagonal();
        }

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs2(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS)
        {
            hs2(V, F, HS);
            auto area = HS.diagonal();
            area = 1.0 / area.array();
            HS = area.asDiagonal();
        }

    }
}
