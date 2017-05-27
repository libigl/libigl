
#ifndef IGL_DISCRETE_EXTERIOR_CALCULUS_H
#define IGL_DISCRETE_EXTERIOR_CALCULUS_H

#include "igl_inline.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
    namespace dec2d
    {
        /*
         * Discrete Exterior Calculus Operators based on Halfedge concept
         *
         * There are two operators in exterior calculus: derivative operator d, and 
         * hodge star operator *. They have different behaviors when operate on 
         * different object (0-form, 1-from and 2-from), so here we distinguash them
         * by adding prefix
         *
         *         V  --- d0 -->  E  -- d1 -->  F       (original graph)
         *         |              |             | 
         *      *0 |           *1 |          *2 |
         *         v              v             v
         *         V* <-- ~d0 --- E* <-- ~d1 -- F*      (dual graph)
         *
         *   V: Represents of 0-form on vertex. 
         *      R^|V|, Row i represents 0-form value on v_i.
         *
         *   E: Represents of 1-form on halfedge.
         *      R^|E| = R^|3F|, Row 3*i+j represents integration of 1-form along 
         *      the j th halfedge of orinted triangle i, in its interior.
         *
         *   F: Represents of 2-form on triangle face.
         *      R^|F|, row i represents integration of 2-form on triangle i.
         *   
         *   V*: Represents of 2-form on dual vertex.
         *      R^|V|, row i represents integration of 2-form on dual area of vertex i.
         *
         *   E*: Represents of 1-form on dual halfedge.
         *      R^|E| = R^|3F|, row 3*i+j represents integration of 1-form along
         *      the dual edge of j th halfedge of orinted triangle i.
         *
         *   F*: Represents of 0-form on dual triangle face.
         *      R^|F|, row i represents 0-form value on dual vertex of triangle i.
         *
         *   d0:  d operator on 0-form represented by vertex, results in 1-form 
         *      represented by halfedge. R^{|E|*|V|}.
         *   d1:  d operator on 1-form represented by halfedge, results in 2-form 
         *      represented by triangle face. R^{|F|*|E|}.
         *   ~d0: d operator on 1-form represented by dual halfedge, results in 
         *      2-form represented by dual vertex. R^{|V|*|E|}.
         *   ~d1: d operator on 0-form represented by dual face, results in 1-form 
         *      represented by dual halfedge. R^{|E|*|F|}.
         *
         *   *0: hodge-star operator on 0-form represented by vertex, results in 
         *      2-form represented by dual vertex, R^{|V|*|V|}.
         *   *1: hodge-star operator on 1-form represented by halfedge, results 
         *      in 1-form represented by dual halfedge, R^{|E|*|E|}.
         *   *2: hodge-star operator on 2-form represented by triangle face, 
         *      results in 0-form represented by dual triangle face, R^{|F|*|F|}.
         *   ~*0: hodge-star operator on 2-form represented by dual vertex, results 
         *      in 0-form represented by vertex, R^{|V|*|V|}. 
         *   ~*1: hodge-star operator on 1-form represented by dual halfedge, results 
         *      in 1-form represented by halfedge, R^{|E|*|E|}.
         *   ~*2: hodge-star operator on 0-form represented by dual face vertex, 
         *      results in 2-form represented by triangle face, R^{|F|*|F|}.
         */


        /*  Consctruct exterior derivative operator
         *
         *  Templates:
         *      DerivedV derived from vertex positions matrix type: i.e. MatrixXd
         *      DerivedF derived from face indices matrix type: i.e. MatrixXi
         *      Scalar derived from sparse matrix type: i.e. SparseMatrix<double>
         *
         *  Inputs:
         *      V eigen matrix #V by 3
         *      F #F by 3 list of mesh faces (must be triangles)
         *
         *  Outputs:
         *      D exterior derivative operator represented in SparseMatrix
         */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void d0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D);
        
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void d1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D);

        /* ~d0 = - 1/2 * transpose(d0) */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_d0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D);

        /* ~d1 = 2 * transpose(d1) */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_d1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& D);
        
        /*  Consctruct hodge star operator
         *
         *  Templates:
         *      DerivedV derived from vertex positions matrix type: i.e. MatrixXd
         *      DerivedF derived from face indices matrix type: i.e. MatrixXi
         *      Scalar derived from sparse matrix type: i.e. SparseMatrix<double>
         *
         *  Inputs:
         *      V eigen matrix #V by 3
         *      F #F by 3 list of mesh faces (must be triangles)
         *
         *  Outputs:
         *      HS hodge star operator represented in SparseMatrix
         */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);

        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void hs2(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);

        /* ~*0 = inverse(*0) */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs0(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);

        /* ~*1 = - inverse(*1) */
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs1(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);
        
        /* ~*2 = inverse(*2)*/
        template <typename DerivedV, typename DerivedF, typename Scalar>
        IGL_INLINE void dual_hs2(
                const Eigen::MatrixBase<DerivedV>& V,
                const Eigen::MatrixBase<DerivedF>& F,
                Eigen::SparseMatrix<Scalar>& HS);
    }
}


#ifndef IGL_STATIC_LIBRARY
#include "discrete_exterior_calculus.cpp"
#endif

#endif
