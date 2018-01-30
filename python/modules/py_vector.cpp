// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include "../python_shared.h"

/// Creates Python bindings for a dynamic Eigen matrix
template <typename Type>
py::class_<Type> bind_eigen_2(py::module &m, const char *name) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> matrix(m, name, py::buffer_protocol());

    matrix
        /* Constructors */
        .def(py::init<>())
        .def(py::init<size_t, size_t>())
        .def("__init__", [](Type &m, Scalar f) {
            new (&m) Type(1, 1);
            m(0, 0) = f;
        })
        .def("__init__", [](Type &m, py::buffer b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<Scalar>::format())
                throw std::runtime_error("Incompatible buffer format!");
            if (info.ndim == 1) {
                new (&m) Type(info.shape[0], 1);
                memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
            } else if (info.ndim == 2) {
                if (info.strides[0] == sizeof(Scalar)) {
                    new (&m) Type(info.shape[0], info.shape[1]);
                    memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
                } else {
                    new (&m) Type(info.shape[1], info.shape[0]);
                    memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
                    m.transposeInPlace();
                }
            } else {
                throw std::runtime_error("Incompatible buffer dimension!");
            }
        })
        .def("__init__", [](Type &m, std::vector<std::vector< Scalar> >& b) {
          if (b.size() == 0)
          {
            new (&m) Type(0, 0);
            return;
          }

          // Size checks
          unsigned rows = b.size();
          unsigned cols = b[0].size();
          for (unsigned i=0;i<rows;++i)
            if (b[i].size() != cols)
              throw std::runtime_error("All rows should have the same size!");

          new (&m) Type(rows, cols);

          m.resize(rows,cols);
          for (unsigned i=0;i<rows;++i)
            for (unsigned j=0;j<cols;++j)
              m(i,j) = b[i][j];

          return;
        })
        .def("__init__", [](Type &m, std::vector<Scalar>& b) {
          if (b.size() == 0)
          {
            new (&m) Type(0, 0);
            return;
          }

          // Size checks
          unsigned rows = b.size();
          unsigned cols = 1;

          new (&m) Type(rows, cols);

          m.resize(rows,cols);
          for (unsigned i=0;i<rows;++i)
            m(i,0) = b[i];

          return;
        })


        /* Size query functions */
        .def("size", [](const Type &m) { return m.size(); })
        .def("cols", [](const Type &m) { return m.cols(); })
        .def("rows", [](const Type &m) { return m.rows(); })
        .def("shape", [](const Type &m) { return std::tuple<int,int>(m.rows(), m.cols()); })

        /* Extract rows and columns */
        .def("col", [](const Type &m, int i) {
            if (i<0 || i>=m.cols())
              throw std::runtime_error("Column index out of bound.");
            return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(m.col(i));
        })
        .def("row", [](const Type &m, int i) {
            if (i<0 || i>=m.rows())
              throw std::runtime_error("Row index out of bound.");
            return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(m.row(i));
        })


        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setIdentity", [](Type &m) { m.setIdentity(); })
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })
        .def("setRandom", [](Type &m) { m.setRandom(); })

        .def("setZero", [](Type &m, const int& r, const int& c) { m.setZero(r,c); })
        .def("setIdentity", [](Type &m, const int& r, const int& c) { m.setIdentity(r,c); })
        .def("setConstant", [](Type &m, const int& r, const int& c, Scalar value) { m.setConstant(r,c,value); })
        .def("setRandom", [](Type &m, const int& r, const int& c) { m.setRandom(r,c); })

        .def("setCol", [](Type &m, int i, const Type& v) { m.col(i) = v; })
        .def("setRow", [](Type &m, int i, const Type& v) { m.row(i) = v; })

        .def("setBlock", [](Type &m, int i, int j, int p, int q, const Type& v) { m.block(i,j,p,q) = v; })
        .def("block", [](Type &m, int i, int j, int p, int q) { return Type(m.block(i,j,p,q)); })

        .def("rightCols", [](Type &m, const int& k) { return Type(m.rightCols(k)); })
        .def("leftCols", [](Type &m, const int& k) { return Type(m.leftCols(k)); })

        .def("setLeftCols", [](Type &m, const int& k, const Type& v) { return Type(m.leftCols(k) = v); })
        .def("setRightCols", [](Type &m, const int& k, const Type& v) { return Type(m.rightCols(k) = v); })

        .def("topRows", [](Type &m, const int& k) { return Type(m.topRows(k)); })
        .def("bottomRows", [](Type &m, const int& k) { return Type(m.bottomRows(k)); })

        .def("setTopRows", [](Type &m, const int& k, const Type& v) { return Type(m.topRows(k) = v); })
        .def("setBottomRows", [](Type &m, const int& k, const Type& v) { return Type(m.bottomRows(k) = v); })

        .def("topLeftCorner", [](Type &m, const int& p, const int&q) { return Type(m.topLeftCorner(p,q)); })
        .def("bottomLeftCorner", [](Type &m, const int& p, const int&q) { return Type(m.bottomLeftCorner(p,q)); })
        .def("topRightCorner", [](Type &m, const int& p, const int&q) { return Type(m.topRightCorner(p,q)); })
        .def("bottomRightCorner", [](Type &m, const int& p, const int&q) { return Type(m.bottomRightCorner(p,q)); })

        /* Resizing */
        .def("resize", [](Type &m, size_t s0, size_t s1) { m.resize(s0, s1); })
        .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        .def("conservativeResize", [](Type &m, size_t s0, size_t s1) { m.conservativeResize(s0, s1); })


        .def("mean", [](const Type &m) {return m.mean();})

        .def("sum", [](const Type &m) {return m.sum();})
        .def("prod", [](const Type &m) {return m.prod();})
        .def("trace", [](const Type &m) {return m.trace();})
        .def("norm", [](const Type &m) {return m.norm();})
        .def("squaredNorm", [](const Type &m) {return m.squaredNorm();})
        .def("squaredMean", [](const Type &m) {return m.array().square().mean();})

        .def("minCoeff", [](const Type &m) {return m.minCoeff();} )
        .def("maxCoeff", [](const Type &m) {return m.maxCoeff();} )

        .def("castdouble", [](const Type &m) {return Eigen::MatrixXd(m.template cast<double>());})
        .def("castint", [](const Type &m) {return Eigen::MatrixXi(m.template cast<int>());})

        /* Component-wise operations */
        .def("cwiseAbs", &Type::cwiseAbs)
        .def("cwiseAbs2", &Type::cwiseAbs2)
        .def("cwiseSqrt", &Type::cwiseSqrt)
        .def("cwiseInverse", &Type::cwiseInverse)
        .def("cwiseMin", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMin(m2); })
        .def("cwiseMax", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMax(m2); })
        .def("cwiseMin", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMin(s); })
        .def("cwiseMax", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMax(s); })
        .def("cwiseProduct", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseProduct(m2); })
        .def("cwiseQuotient", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseQuotient(m2); })

        /* Row and column-wise operations */
        .def("rowwiseSet", [](Type &m, const Type &m2) {return Type(m.rowwise() = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>(m2));} )
        .def("rowwiseSum", [](const Type &m) {return Type(m.rowwise().sum());} )
        .def("rowwiseProd", [](const Type &m) {return Type(m.rowwise().prod());} )
        .def("rowwiseMean", [](const Type &m) {return Type(m.rowwise().mean());} )
        .def("rowwiseNorm", [](const Type &m) {return Type(m.rowwise().norm());} )
        .def("rowwiseNormalized", [](const Type &m) {return Type(m.rowwise().normalized());} )
        .def("rowwiseReverse", [](const Type &m) {return Type(m.rowwise().reverse());} )
        .def("rowwiseMinCoeff", [](const Type &m) {return Type(m.rowwise().minCoeff());} )
        .def("rowwiseMaxCoeff", [](const Type &m) {return Type(m.rowwise().maxCoeff());} )

        .def("colwiseSet", [](Type &m, const Type &m2) {return Type(m.colwise() = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>(m2));} )
        .def("colwiseSum", [](const Type &m) {return Type(m.colwise().sum());} )
        .def("colwiseProd", [](const Type &m) {return Type(m.colwise().prod());} )
        .def("colwiseMean", [](const Type &m) {return Type(m.colwise().mean());} )
        .def("colwiseNorm", [](const Type &m) {return Type(m.colwise().norm());} )
        .def("colwiseNormalized", [](const Type &m) {return Type(m.colwise().normalized());} )
        .def("colwiseReverse", [](const Type &m) {return Type(m.colwise().reverse());} )
        .def("colwiseMinCoeff", [](const Type &m) {return Type(m.colwise().minCoeff());} )
        .def("colwiseMaxCoeff", [](const Type &m) {return Type(m.colwise().maxCoeff());} )

        .def("replicate", [](const Type &m, const int& r, const int& c) {return Type(m.replicate(r,c));} )
        .def("asDiagonal", [](const Type &m) {return Eigen::DiagonalMatrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(m.asDiagonal());} )

        .def("sparseView", [](Type &m) { return Eigen::SparseMatrix<Scalar>(m.sparseView()); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * py::self)
        // .def_cast(py::self - Scalar())
        // .def_cast(py::self * Scalar())
        // .def_cast(py::self / Scalar())

        .def("__mul__", []
        (const Type &a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a * b);
        })
        .def("__rmul__", [](const Type& a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b * a);
        })

        .def("__add__", []
        (const Type &a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a.array() + b);
        })
        .def("__radd__", [](const Type& a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b + a.array());
        })

        .def("__sub__", []
        (const Type &a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a.array() - b);
        })
        .def("__rsub__", [](const Type& a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b - a.array());
        })

        .def("__div__", []
        (const Type &a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a / b);
        })

        .def("__truediv__", []
        (const Type &a, const Scalar& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a / b);
        })

        /* Arithmetic in-place operators */
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__lt__", []
        (const Type &a, const Scalar& b) -> Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>
        {
          return Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>(a.array() < b);
        })
        .def("__gt__", []
        (const Type &a, const Scalar& b) -> Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>
        {
          return Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>(a.array() > b);
        })
        .def("__le__", []
        (const Type &a, const Scalar& b) -> Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>
        {
          return Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>(a.array() <= b);
        })
        .def("__ge__", []
        (const Type &a, const Scalar& b) -> Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>
        {
          return Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>(a.array() >= b);
        })

        .def("transposeInPlace", [](Type &m) { m.transposeInPlace(); })
        /* Other transformations */
        .def("transpose", [](Type &m) -> Type { return m.transpose(); })
        /* Python protocol implementations */
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        })
        .def("__getitem__", [](const Type &m, std::pair<size_t, size_t> i) {
            if (i.first >= (size_t) m.rows() || i.second >= (size_t) m.cols())
                throw py::index_error();
            return m(i.first, i.second);
         })
        .def("__setitem__", [](Type &m, std::pair<size_t, size_t> i, Scalar v) {
            if (i.first >= (size_t) m.rows() || i.second >= (size_t) m.cols())
                throw py::index_error();
            m(i.first, i.second) = v;
         })

         .def("__getitem__", [](const Type &m, size_t i) {
             if (i >= (size_t) m.size())
                 throw py::index_error();
             return m(i);
          })
         .def("__setitem__", [](Type &m, size_t i, Scalar v) {
           if (i >= (size_t) m.size())
                 throw py::index_error();
             m(i) = v;
          })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                /* Pointer to buffer */
                sizeof(Scalar),          /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::format(),
                2,                       /* Number of dimensions */
                { (size_t) m.rows(),     /* Buffer dimensions */
                  (size_t) m.cols() },
                { sizeof(Scalar),        /* Strides (in bytes) for each index */
                  sizeof(Scalar) * m.rows() }
            );
         })

        /* Static initializers */
        .def_static("Zero", [](size_t n, size_t m) { return Type(Type::Zero(n, m)); })
        .def_static("Random", [](size_t n, size_t m) { return Type(Type::Random(n, m)); })
        .def_static("Ones", [](size_t n, size_t m) { return Type(Type::Ones(n, m)); })
        .def_static("Constant", [](size_t n, size_t m, Scalar value) { return Type(Type::Constant(n, m, value)); })
        .def_static("Identity", [](size_t n, size_t m) { return Type(Type::Identity(n, m)); })
        .def("MapMatrix", [](const Type& m, size_t r, size_t c)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(m.data(),r,c));
        })

        .def("copy", [](const Type &m) { return Type(m); })

        ;
    return matrix;
}

/// Creates Python bindings for a dynamic Eigen sparse order-2 tensor (i.e. a matrix)
template <typename Type>
py::class_<Type> bind_eigen_sparse_2(py::module &m, const char *name) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> matrix(m, name, py::buffer_protocol());

    matrix
        /* Constructors */
        .def(py::init<>())
        .def(py::init<size_t, size_t>())
        // .def("__init__", [](Type &m, Scalar f) {
        //     new (&m) Type(1, 1);
        //     m(0, 0) = f;
        // })
        // .def("__init__", [](Type &m, py::buffer b) {
        //     py::buffer_info info = b.request();
        //     if (info.format != py::format_descriptor<Scalar>::value())
        //         throw std::runtime_error("Incompatible buffer format!");
        //     if (info.ndim == 1) {
        //         new (&m) Type(info.shape[0], 1);
        //         memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
        //     } else if (info.ndim == 2) {
        //         if (info.strides[0] == sizeof(Scalar)) {
        //             new (&m) Type(info.shape[0], info.shape[1]);
        //             memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
        //         } else {
        //             new (&m) Type(info.shape[1], info.shape[0]);
        //             memcpy(m.data(), info.ptr, sizeof(Scalar) * m.size());
        //             m.transposeInPlace();
        //         }
        //     } else {
        //         throw std::runtime_error("Incompatible buffer dimension!");
        //     }
        // })

        /* Size query functions */
        .def("size", [](const Type &m) { return m.size(); })
        .def("cols", [](const Type &m) { return m.cols(); })
        .def("rows", [](const Type &m) { return m.rows(); })
        .def("shape", [](const Type &m) { return std::tuple<int,int>(m.rows(), m.cols()); })


        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setIdentity", [](Type &m) { m.setIdentity(); })

        .def("transpose", [](Type &m) { return Type(m.transpose()); })
        .def("norm", [](Type &m) { return m.norm(); })

        /* Resizing */
        // .def("resize", [](Type &m, size_t s0, size_t s1) { m.resize(s0, s1); })
        // .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        // .def("conservativeResize", [](Type &m, size_t s0, size_t s1) { m.conservativeResize(s0, s1); })

        /* Component-wise operations */
        // .def("cwiseAbs", &Type::cwiseAbs)
        // .def("cwiseAbs2", &Type::cwiseAbs2)
        // .def("cwiseSqrt", &Type::cwiseSqrt)
        // .def("cwiseInverse", &Type::cwiseInverse)
        // .def("cwiseMin", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMin(m2); })
        // .def("cwiseMax", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseMax(m2); })
        // .def("cwiseMin", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMin(s); })
        // .def("cwiseMax", [](const Type &m1, Scalar s) -> Type { return m1.cwiseMax(s); })
        // .def("cwiseProduct", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseProduct(m2); })
        // .def("cwiseQuotient", [](const Type &m1, const Type &m2) -> Type { return m1.cwiseQuotient(m2); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * py::self)
        .def_cast(py::self * Scalar())
        .def_cast(Scalar() * py::self)
        // Special case, sparse * dense produces a dense matrix

        // .def("__mul__", []
        // (const Type &a, const Scalar& b)
        // {
        //   return Type(a * b);
        // })
        // .def("__rmul__", [](const Type& a, const Scalar& b)
        // {
        //   return Type(b * a);
        // })

        .def("__mul__", []
        (const Type &a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a * b);
        })
        .def("__rmul__", [](const Type& a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b * a);
        })

        .def("__mul__", []
        (const Type &a, const Eigen::DiagonalMatrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Type(a * b);
        })
        .def("__rmul__", [](const Type& a, const Eigen::DiagonalMatrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Type(b * a);
        })

        //.def(py::self * Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>())
//        .def_cast(py::self / Scalar())

        /* Arithmetic in-place operators */
        // .def_cast(py::self += py::self)
        // .def_cast(py::self -= py::self)
        // .def_cast(py::self *= py::self)
        // .def_cast(py::self *= Scalar())
        // .def_cast(py::self /= Scalar())

        /* Comparison operators */
        // .def(py::self == py::self)
        // .def(py::self != py::self)

        // .def("transposeInPlace", [](Type &m) { m.transposeInPlace(); })
        // /* Other transformations */
        // .def("transpose", [](Type &m) -> Type { return m.transpose(); })

        /* Python protocol implementations */
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        })

        /* Static initializers */
        // .def_static("Zero", [](size_t n, size_t m) { return Type(Type::Zero(n, m)); })
        // .def_static("Ones", [](size_t n, size_t m) { return Type(Type::Ones(n, m)); })
        // .def_static("Constant", [](size_t n, size_t m, Scalar value) { return Type(Type::Constant(n, m, value)); })
        // .def_static("Identity", [](size_t n, size_t m) { return Type(Type::Identity(n, m)); })
        .def("toCOO",[](const Type& m)
        {
          Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> t(m.nonZeros(),3);
          int count = 0;
          for (int k=0; k<m.outerSize(); ++k)
            for (typename Type::InnerIterator it(m,k); it; ++it)
              t.row(count++) << it.row(), it.col(), it.value();
          return t;
        })

        .def("fromCOO",[](Type& m, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& t, int rows, int cols)
        {
          typedef Eigen::Triplet<Scalar> T;
          std::vector<T> tripletList;
          tripletList.reserve(t.rows());
          for(unsigned i=0;i<t.rows();++i)
            tripletList.push_back(T(round(t(i,0)),round(t(i,1)),t(i,2)));

          if (rows == -1)
            rows = t.col(0).maxCoeff()+1;

          if (cols == -1)
            cols = t.col(1).maxCoeff()+1;

          m.resize(rows,cols);
          m.setFromTriplets(tripletList.begin(), tripletList.end());
        }, py::arg("t"), py::arg("rows") = -1, py::arg("cols") = -1)

        .def("insert",[](Type& m, const int row, const int col, const Scalar value)
        {
          return m.insert(row,col) = value;
        }, py::arg("row"), py::arg("col"), py::arg("value"))

        .def("makeCompressed",[](Type& m)
        {
          return m.makeCompressed();
        })

        .def("diagonal", [](const Type &m) {return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(m.diagonal());} )

        ;
    return matrix;
}

/// Creates Python bindings for a diagonal Eigen sparse order-2 tensor (i.e. a matrix)
template <typename Type>
py::class_<Type> bind_eigen_diagonal_2(py::module &m, const char *name) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> matrix(m, name, py::buffer_protocol());

    matrix
        /* Constructors */
        .def(py::init<>())
        //.def(py::init<size_t, size_t>())

        /* Size query functions */
        .def("size", [](const Type &m) { return m.size(); })
        .def("cols", [](const Type &m) { return m.cols(); })
        .def("rows", [](const Type &m) { return m.rows(); })
        .def("shape", [](const Type &m) { return std::tuple<int,int>(m.rows(), m.cols()); })

        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setIdentity", [](Type &m) { m.setIdentity(); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        // .def_cast(-py::self)
        // .def_cast(py::self + py::self)
        // .def_cast(py::self - py::self)
        // .def_cast(py::self * py::self)
        .def_cast(py::self * Scalar())
        .def_cast(Scalar() * py::self)

        // // Special case, sparse * dense produces a dense matrix
        // .def("__mul__", []
        // (const Type &a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        // {
        //   return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a * b);
        // })
        // .def("__rmul__", [](const Type& a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        // {
        //   return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b * a);
        // })

        .def("__mul__", []
        (const Type &a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a * b);
        })
        .def("__rmul__", [](const Type& a, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& b)
        {
          return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b * a);
        })

        .def("__mul__", []
        (const Type &a, const Eigen::SparseMatrix<Scalar>& b)
        {
          return Eigen::SparseMatrix<Scalar>(a * b);
        })
        .def("__rmul__", [](const Type& a, const Eigen::SparseMatrix<Scalar>& b)
        {
          return Eigen::SparseMatrix<Scalar>(b * a);
        })

        /* Python protocol implementations */
        .def("__repr__", [](const Type &/*v*/) {
            std::ostringstream oss;
            oss << "<< operator undefined for diagonal matrices";
            return oss.str();
        })

        /* Other transformations */

        ;
    return matrix;
}


void python_export_vector(py::module &m) {

  py::module me = m.def_submodule(
    "eigen", "Wrappers for Eigen types");

    /* Bindings for VectorXd */
    // bind_eigen_1<Eigen::VectorXd> (me, "VectorXd");
    // py::implicitly_convertible<py::buffer, Eigen::VectorXd>();
    // py::implicitly_convertible<double, Eigen::VectorXd>();

    /* Bindings for VectorXi */
    // bind_eigen_1<Eigen::VectorXi> (me, "VectorXi");
    // py::implicitly_convertible<py::buffer, Eigen::VectorXi>();
    // py::implicitly_convertible<double, Eigen::VectorXi>();

    /* Bindings for MatrixXd */
    bind_eigen_2<Eigen::MatrixXd> (me, "MatrixXd");
    //py::implicitly_convertible<py::buffer, Eigen::MatrixXd>();
    //py::implicitly_convertible<double, Eigen::MatrixXd>();

    /* Bindings for MatrixXi */
    bind_eigen_2<Eigen::MatrixXi> (me, "MatrixXi");
    //py::implicitly_convertible<py::buffer, Eigen::MatrixXi>();
    //py::implicitly_convertible<double, Eigen::MatrixXi>();

    /* Bindings for MatrixXb */
    bind_eigen_2<Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> > (me, "MatrixXb");

    /* Bindings for MatrixXuc */
    bind_eigen_2<Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> > (me, "MatrixXuc");
    // py::implicitly_convertible<py::buffer, Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> >();
    // py::implicitly_convertible<double, Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> >();

    // /* Bindings for Vector3d */
    // auto vector3 = bind_eigen_1_3<Eigen::Vector3d>(me, "Vector3d");
    // vector3
    //     .def("norm", [](const Eigen::Vector3d &v) { return v.norm(); })
    //     .def("squaredNorm", [](const Eigen::Vector3d &v) { return v.squaredNorm(); })
    //     .def("normalize", [](Eigen::Vector3d &v) { v.normalize(); })
    //     .def("normalized", [](const Eigen::Vector3d &v) -> Eigen::Vector3d { return v.normalized(); })
    //     .def("dot", [](const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) { return v1.dot(v2); })
    //     .def("cross", [](const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) -> Eigen::Vector3d { return v1.cross(v2); })
    //     .def_property("x", [](const Eigen::Vector3d &v) -> double { return v.x(); },
    //                        [](Eigen::Vector3d &v, double x) { v.x() = x; }, "X coordinate")
    //     .def_property("y", [](const Eigen::Vector3d &v) -> double { return v.y(); },
    //                        [](Eigen::Vector3d &v, double y) { v.y() = y; }, "Y coordinate")
    //     .def_property("z", [](const Eigen::Vector3d &v) -> double { return v.z(); },
    //                        [](Eigen::Vector3d &v, double z) { v.z() = z; }, "Z coordinate");
    //
    // py::implicitly_convertible<py::buffer, Eigen::Vector3d>();
    // py::implicitly_convertible<double, Eigen::Vector3d>();

    /* Bindings for SparseMatrix<double> */
    bind_eigen_sparse_2< Eigen::SparseMatrix<double> > (me, "SparseMatrixd");

    /* Bindings for SparseMatrix<int> */
    bind_eigen_sparse_2< Eigen::SparseMatrix<int> > (me, "SparseMatrixi");

    /* Bindings for DiagonalMatrix<double> */
    bind_eigen_diagonal_2< Eigen::DiagonalMatrix<double,Eigen::Dynamic,Eigen::Dynamic> > (me, "DiagonalMatrixd");

    /* Bindings for DiagonalMatrix<int> */
    bind_eigen_diagonal_2< Eigen::DiagonalMatrix<int,Eigen::Dynamic,Eigen::Dynamic> > (me, "DiagonalMatrixi");

    /* Bindings for SimplicialLLT*/
    py::class_<Eigen::SimplicialLLT<Eigen::SparseMatrix<double > >> simpliciallltsparse(me, "SimplicialLLTsparse");

    simpliciallltsparse
    .def(py::init<>())
    .def(py::init<Eigen::SparseMatrix<double>>())
    .def("info",[](const Eigen::SimplicialLLT<Eigen::SparseMatrix<double > >& s)
    {
       if (s.info() == Eigen::Success)
          return "Success";
       else
          return "Numerical Issue";
    })
    .def("analyzePattern",[](Eigen::SimplicialLLT<Eigen::SparseMatrix<double > >& s, const Eigen::SparseMatrix<double>& a) { return s.analyzePattern(a); })
    .def("factorize",[](Eigen::SimplicialLLT<Eigen::SparseMatrix<double > >& s, const Eigen::SparseMatrix<double>& a) { return s.factorize(a); })
    .def("solve",[](const Eigen::SimplicialLLT<Eigen::SparseMatrix<double > >& s, const Eigen::MatrixXd& rhs) { return Eigen::MatrixXd(s.solve(rhs)); })
    ;

    // Bindings for Affine3d
    py::class_<Eigen::Affine3d > affine3d(me, "Affine3d");

    affine3d
    .def(py::init<>())
    .def_static("Identity", []() { return Eigen::Affine3d::Identity(); })
    .def("setIdentity",[](Eigen::Affine3d& a){
        return a.setIdentity();
    })
    .def("rotate",[](Eigen::Affine3d& a, double angle, Eigen::MatrixXd axis) {
        assert_is_Vector3("axis", axis);
        return a.rotate(Eigen::AngleAxisd(angle, Eigen::Vector3d(axis)));
    })
    .def("rotate",[](Eigen::Affine3d& a, Eigen::Quaterniond quat) {
        return a.rotate(quat);
    })
    .def("translate",[](Eigen::Affine3d& a, Eigen::MatrixXd offset) {
        assert_is_Vector3("offset", offset);
        return a.translate(Eigen::Vector3d(offset));
    })
    .def("matrix", [](Eigen::Affine3d& a) -> Eigen::MatrixXd {
        return Eigen::MatrixXd(a.matrix());
    })
    ;
    // Bindings for Quaterniond
    py::class_<Eigen::Quaterniond > quaterniond(me, "Quaterniond");
    
    quaterniond
    .def(py::init<>())
    .def(py::init<double, double, double, double>())
    .def("__init__", [](Eigen::Quaterniond &q, double angle, Eigen::MatrixXd axis) {
        assert_is_Vector3("axis", axis);
        new (&q) Eigen::Quaterniond(Eigen::AngleAxisd(angle, Eigen::Vector3d(axis)));
    })
    .def_static("Identity", []() { return Eigen::Quaterniond::Identity(); })
    .def("__repr__", [](const Eigen::Quaterniond &v) {
        std::ostringstream oss;
        oss << "(" << v.w() << ", " << v.x() << ", " << v.y() << ", " << v.z() << ")";
        return oss.str();
    })
    .def("conjugate",[](Eigen::Quaterniond& q) {
        return q.conjugate();
    })
    .def("normalize",[](Eigen::Quaterniond& q) {
        return q.normalize();
    })
    .def("slerp",[](Eigen::Quaterniond& q, double & t, Eigen::Quaterniond other) {
        return q.slerp(t, other);
    })
//    .def_cast(-py::self)
//    .def_cast(py::self + py::self)
//    .def_cast(py::self - py::self)
    .def_cast(py::self * py::self)
    // .def_cast(py::self - Scalar())
    // .def_cast(py::self * Scalar())
    // .def_cast(py::self / Scalar())

//    .def("__mul__", []
//    (const Type &a, const Scalar& b)
//    {
//      return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(a * b);
//    })
//    .def("__rmul__", [](const Type& a, const Scalar& b)
//    {
//      return Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(b * a);
//    })
    ;

    




}
