#include <Eigen/Dense>

#include "python.h"

template <typename Type> void init_fixed_from_buffer_3(Type &v, py::buffer &b) {
    typedef typename Type::Scalar Scalar;

    py::buffer_info info = b.request();
    if (info.format != py::format_descriptor<Scalar>::value())
        throw std::runtime_error("Incompatible buffer format!");
    if (!((info.ndim == 1 && info.strides[0] == sizeof(Scalar)) ||
          (info.ndim == 2 &&
              ((info.shape[0] == 1 && info.strides[0] == sizeof(Scalar) &&
                info.shape[1] == 3) ||
               (info.shape[1] == 1 && info.strides[1] == sizeof(Scalar) &&
                info.shape[0] == 3)))))
        throw std::runtime_error("Incompatible buffer dimension!");

    memcpy(v.data(), info.ptr, sizeof(Scalar) * 3);
}

/// Creates Python bindings for an Eigen order-1 tensor of size 3 (i.e. a vector/normal/point)
template <typename Type>
py::class_<Type> bind_eigen_1_3(py::module &m, const char *name,
                                py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    py::class_<Type> vector(m, name, parent);
    vector
        /* Constructors */
        .def(py::init<>())
        .def(py::init<Scalar>())
        .def(py::init<Scalar, Scalar, Scalar>())
        .def("__init__", [](Type &v, const std::vector<Scalar> &v2) {
            if (v2.size() != 3)
                throw std::runtime_error("Incompatible size!");
            memcpy(v.data(), &v2[0], sizeof(Scalar) * 3);
        })
        .def("__init__", [](Type &v, py::buffer b) {
            init_fixed_from_buffer_3(v, b);
        })

        /* Initialization */
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })
        .def("setZero", [](Type &m) { m.setZero(); })

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Matrix to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

        /* Python protocol implementations */
        .def("__len__", [](const Type &) { return (int) 3; })
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v;
            return oss.str();
        })
        .def("__getitem__", [](const Type &c, int i) {
            if (i < 0 || i >= 3)
                throw py::index_error();
            return c[i];
         })
        .def("__setitem__", [](Type &c, int i, Scalar v) {
             if (i < 0 || i >= 3)
                 throw py::index_error();
            c[i] = v;
         })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),        /* Pointer to buffer */
                sizeof(Scalar),  /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                1, { (size_t) 3 },
                { sizeof(Scalar) }
            );
        });
    return vector;
}

/// Creates Python bindings for a dynamic Eigen order-1 tensor (i.e. a vector)
template <typename Type>
py::class_<Type> bind_eigen_1(py::module &m, const char *name,
                              py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> vector(m, name, parent);

    vector
        /* Constructors */
        .def(py::init<>())
        .def(py::init<size_t>())
        .def("__init__", [](Type &v, const std::vector<Scalar> &v2) {
            new (&v) Type(v2.size());
            memcpy(v.data(), &v2[0], sizeof(Scalar) * v2.size());
        })
        .def("__init__", [](Type &v, py::buffer b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<Scalar>::value()) {
                throw std::runtime_error("Incompatible buffer format!");
            } else if (info.ndim == 1 && info.strides[0] == sizeof(Scalar)) {
                new (&v) Type(info.shape[0]);
                memcpy(v.data(), info.ptr, sizeof(Scalar) * info.shape[0]);
            } else if (info.ndim == 2 && ((info.shape[0] == 1 && info.strides[0] == sizeof(Scalar))
                                       || (info.shape[1] == 1 && info.strides[1] == sizeof(Scalar)))) {
                new (&v) Type(info.shape[0] * info.shape[1]);
                memcpy(v.data(), info.ptr, sizeof(Scalar) * info.shape[0] * info.shape[1]);
            } else {
                throw std::runtime_error("Incompatible buffer dimension!");
            }
        })

        /* Size query functions */
        .def("size", [](const Type &m) { return m.size(); })
        .def("cols", &Type::cols)
        .def("rows", &Type::rows)

        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })

        /* Resizing */
        .def("resize", [](Type &m, size_t s0) { m.resize(s0); })
        .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        .def("conservativeResize", [](Type &m, size_t s0) { m.conservativeResize(s0); })

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

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())

        /* Arithmetic in-place operators */
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

        /* Python protocol implementations */
        .def("__repr__", [](const Type &v) {
            std::ostringstream oss;
            oss << v.transpose();
            return oss.str();
        })
        .def("__getitem__", [](const Type &m, size_t i) {
            if (i >= (size_t) m.size())
                throw py::index_error();
            return m[i];
         })
        .def("__setitem__", [](Type &m, size_t i, Scalar v) {
            if (i >= (size_t) m.size())
                throw py::index_error();
            m[i] = v;
         })

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                /* Pointer to buffer */
                sizeof(Scalar),          /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                1,                       /* Number of dimensions */
                { (size_t) m.size() },   /* Buffer dimensions */
                { sizeof(Scalar) }       /* Strides (in bytes) for each index */
            );
         })

        /* Static initializers */
        .def_static("Zero", [](size_t n) { return Type(Type::Zero(n)); })
        .def_static("Ones", [](size_t n) { return Type(Type::Ones(n)); })
        .def_static("Constant", [](size_t n, Scalar value) { return Type(Type::Constant(n, value)); });
    return vector;
}

/// Creates Python bindings for a dynamic Eigen order-2 tensor (i.e. a matrix)
template <typename Type>
py::class_<Type> bind_eigen_2(py::module &m, const char *name,
                                py::object parent = py::object()) {
    typedef typename Type::Scalar Scalar;

    /* Many Eigen functions are templated and can't easily be referenced using
       a function pointer, thus a big portion of the binding code below
       instantiates Eigen code using small anonymous wrapper functions */
    py::class_<Type> matrix(m, name, parent);

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
            if (info.format != py::format_descriptor<Scalar>::value())
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

        /* Size query functions */
        .def("size", [](const Type &m) { return m.size(); })
        .def("cols", &Type::cols)
        .def("rows", &Type::rows)

        /* Initialization */
        .def("setZero", [](Type &m) { m.setZero(); })
        .def("setIdentity", [](Type &m) { m.setIdentity(); })
        .def("setConstant", [](Type &m, Scalar value) { m.setConstant(value); })

        /* Resizing */
        .def("resize", [](Type &m, size_t s0, size_t s1) { m.resize(s0, s1); })
        .def("resizeLike", [](Type &m, const Type &m2) { m.resizeLike(m2); })
        .def("conservativeResize", [](Type &m, size_t s0, size_t s1) { m.conservativeResize(s0, s1); })

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

        /* Arithmetic operators (def_cast forcefully casts the result back to a
           Type to avoid type issues with Eigen's crazy expression templates) */
        .def_cast(-py::self)
        .def_cast(py::self + py::self)
        .def_cast(py::self - py::self)
        .def_cast(py::self * py::self)
        .def_cast(py::self * Scalar())
        .def_cast(py::self / Scalar())

        /* Arithmetic in-place operators */
        .def_cast(py::self += py::self)
        .def_cast(py::self -= py::self)
        .def_cast(py::self *= py::self)
        .def_cast(py::self *= Scalar())
        .def_cast(py::self /= Scalar())

        /* Comparison operators */
        .def(py::self == py::self)
        .def(py::self != py::self)

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

        /* Buffer access for interacting with NumPy */
        .def_buffer([](Type &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                /* Pointer to buffer */
                sizeof(Scalar),          /* Size of one scalar */
                /* Python struct-style format descriptor */
                py::format_descriptor<Scalar>::value(),
                2,                       /* Number of dimensions */
                { (size_t) m.rows(),     /* Buffer dimensions */
                  (size_t) m.cols() },
                { sizeof(Scalar),        /* Strides (in bytes) for each index */
                  sizeof(Scalar) * m.rows() }
            );
         })

        /* Static initializers */
        .def_static("Zero", [](size_t n, size_t m) { return Type(Type::Zero(n, m)); })
        .def_static("Ones", [](size_t n, size_t m) { return Type(Type::Ones(n, m)); })
        .def_static("Constant", [](size_t n, size_t m, Scalar value) { return Type(Type::Constant(n, m, value)); })
        .def_static("Identity", [](size_t n, size_t m) { return Type(Type::Identity(n, m)); });
    return matrix;
}

void python_export_vector(py::module &m) {

  py::module me = m.def_submodule(
    "eigen", "Wrappers for Eigen types");

    bind_eigen_1<Eigen::VectorXd> (me, "VectorXd");
    bind_eigen_1<Eigen::VectorXi> (me, "VectorXi");
    bind_eigen_2<Eigen::MatrixXd> (me, "MatrixXd");
    bind_eigen_2<Eigen::MatrixXi> (me, "MatrixXi");


    /* Bindings for <vector.h> */
    auto vector3 = bind_eigen_1_3<Eigen::Vector3d>(me, "Vector3d");
    vector3
        .def("norm", [](const Eigen::Vector3d &v) { return v.norm(); })
        .def("squaredNorm", [](const Eigen::Vector3d &v) { return v.squaredNorm(); })
        .def("normalize", [](Eigen::Vector3d &v) { v.normalize(); })
        .def("normalized", [](const Eigen::Vector3d &v) -> Eigen::Vector3d { return v.normalized(); })
        .def("dot", [](const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) { return v1.dot(v2); })
        .def("cross", [](const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) -> Eigen::Vector3d { return v1.cross(v2); })
        .def_property("x", [](const Eigen::Vector3d &v) -> double { return v.x(); },
                           [](Eigen::Vector3d &v, double x) { v.x() = x; }, "X coordinate")
        .def_property("y", [](const Eigen::Vector3d &v) -> double { return v.y(); },
                           [](Eigen::Vector3d &v, double y) { v.y() = y; }, "Y coordinate")
        .def_property("z", [](const Eigen::Vector3d &v) -> double { return v.z(); },
                           [](Eigen::Vector3d &v, double z) { v.z() = z; }, "Z coordinate");

    py::implicitly_convertible<py::buffer, Eigen::VectorXd>();
    py::implicitly_convertible<py::buffer, Eigen::MatrixXd>();
    py::implicitly_convertible<py::buffer, Eigen::VectorXi>();
    py::implicitly_convertible<py::buffer, Eigen::MatrixXi>();
    py::implicitly_convertible<py::buffer, Eigen::Vector3d>();

    py::implicitly_convertible<double, Eigen::VectorXd>();
    py::implicitly_convertible<double, Eigen::MatrixXd>();
    py::implicitly_convertible<double, Eigen::VectorXi>();
    py::implicitly_convertible<double, Eigen::MatrixXi>();
    py::implicitly_convertible<double, Eigen::Vector3d>();
}
