#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

typedef std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;
PYBIND11_MAKE_OPAQUE(RotationList)

//typedef std::vector<Eigen::Vector3d> TranslationList;
//PYBIND11_MAKE_OPAQUE(TranslationList);

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>)

