#include <igl/polyvector_field_poisson_reconstruction.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/sparse.h>
#include <igl/repdiag.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/colon.h>

#include <Eigen/Sparse>

template <typename DerivedV, typename DerivedF, typename DerivedSF, typename DerivedS, typename DerivedE>
IGL_INLINE double igl::polyvector_field_poisson_reconstruction(
  const Eigen::PlainObjectBase<DerivedV> &Vcut,
  const Eigen::PlainObjectBase<DerivedF> &Fcut,
  const Eigen::PlainObjectBase<DerivedS> &sol3D_combed,
  Eigen::PlainObjectBase<DerivedSF> &scalars,
  Eigen::PlainObjectBase<DerivedS> &sol3D_recon,
  Eigen::PlainObjectBase<DerivedE> &max_error )
  {
    Eigen::SparseMatrix<typename DerivedV::Scalar> gradMatrix;
    igl::grad(Vcut, Fcut, gradMatrix);

    Eigen::VectorXd FAreas;
    igl::doublearea(Vcut, Fcut, FAreas);
    FAreas = FAreas.array() * .5;

    int nf = FAreas.rows();
    Eigen::SparseMatrix<typename DerivedV::Scalar> M,M1;
    Eigen::VectorXi II = igl::colon<int>(0, nf-1);

    igl::sparse(II, II, FAreas, M1);
    igl::repdiag(M1, 3, M) ;

    int half_degree = sol3D_combed.cols()/3;

    sol3D_recon.setZero(sol3D_combed.rows(),sol3D_combed.cols());

    int numF = Fcut.rows();
    scalars.setZero(Vcut.rows(),half_degree);

    Eigen::SparseMatrix<typename DerivedV::Scalar> Q = gradMatrix.transpose()* M *gradMatrix;

    //fix one point at Ik=fix, value at fixed xk=0
    int fix = 0;
    Eigen::VectorXi Ik(1);Ik<<fix;
    Eigen::VectorXd xk(1);xk<<0;

    //unknown indices
    Eigen::VectorXi Iu(Vcut.rows()-1,1);
    Iu<<igl::colon<int>(0, fix-1),  igl::colon<int>(fix+1,Vcut.rows()-1);

    Eigen::SparseMatrix<typename DerivedV::Scalar> Quu, Quk;
    igl::slice(Q, Iu, Iu, Quu);
    igl::slice(Q, Iu, Ik, Quk);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<typename DerivedV::Scalar> > solver;
    solver.compute(Quu);


    Eigen::VectorXd vec; vec.setZero(3*numF,1);
    for (int i =0; i<half_degree; ++i)
    {
      vec<<sol3D_combed.col(i*3+0),sol3D_combed.col(i*3+1),sol3D_combed.col(i*3+2);
      Eigen::VectorXd b = gradMatrix.transpose()* M * vec;
      Eigen::VectorXd bu = igl::slice(b, Iu);

      Eigen::VectorXd rhs = bu-Quk*xk;
      Eigen::MatrixXd yu = solver.solve(rhs);

      Eigen::VectorXi index = i*Eigen::VectorXi::Ones(Iu.rows(),1);
      igl::slice_into(yu, Iu, index, scalars);scalars(Ik[0],i)=xk[0];
    }

    //    evaluate gradient of found scalar function
    for (int i =0; i<half_degree; ++i)
    {
      Eigen::VectorXd vec_poisson = gradMatrix*scalars.col(i);
      sol3D_recon.col(i*3+0) = vec_poisson.segment(0*numF, numF);
      sol3D_recon.col(i*3+1) = vec_poisson.segment(1*numF, numF);
      sol3D_recon.col(i*3+2) = vec_poisson.segment(2*numF, numF);
    }

    max_error.setZero(numF,1);
    for (int i =0; i<half_degree; ++i)
    {
      Eigen::VectorXd diff = (sol3D_recon.block(0, i*3, numF, 3)-sol3D_combed.block(0, i*3, numF, 3)).rowwise().norm();
      diff = diff.array() / sol3D_combed.block(0, i*3, numF, 3).rowwise().norm().array();
      max_error = max_error.cwiseMax(diff.cast<typename DerivedE::Scalar>());
    }

    return max_error.mean();
  }

  #ifdef IGL_STATIC_LIBRARY
  // Explicit template specialization
  template double igl::polyvector_field_poisson_reconstruction<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 1, 0, -1, 1> >&);
  #endif
