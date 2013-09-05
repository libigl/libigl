#include "matlab_format.h"
#include "STR.h"
#include "find.h"

template <typename DerivedM>
IGL_INLINE const Eigen::WithFormat< DerivedM > igl::matlab_format(
  const Eigen::PlainObjectBase<DerivedM> & M,
  const std::string name = "")
{
  using namespace igl;
  using namespace std;
  string prefix = "";
  if(!name.empty())
  {
    prefix = name + " = ";
  }

  return M.format(Eigen::IOFormat(
    Eigen::FullPrecision,
    0,
    " ",
    "\n",
    "",
    "",
    // This seems like a bit of a hack since I would expect the rows to align
    // with out this extra spacing on the first line
    prefix + "[\n  ",
    "\n];"));
}

template <typename DerivedS>
IGL_INLINE const std::string
igl::matlab_format(
  const Eigen::SparseMatrix<DerivedS> & S,
  const std::string name = "")
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  Matrix<typename Eigen::SparseMatrix<DerivedS>::Scalar,Dynamic,1> I,J,V;
  Matrix<DerivedS,Dynamic,Dynamic> SIJV;
  find(S,I,J,V);
  I.array() += 1;
  J.array() += 1;
  SIJV.resize(V.rows(),3);
  SIJV << I,J,V;
  string prefix = "";
  string suffix = "";
  if(!name.empty())
  {
    prefix = name + "IJV = ";
    suffix = "\n"+name + " = sparse("+name+"IJV(:,1),"+name+"IJV(:,2),"+name+"IJV(:,3));";
  }
  return STR(""<<
    SIJV.format(Eigen::IOFormat(
    Eigen::FullPrecision,
    0,
    " ",
    "\n",
    "",
    "",
    // This seems like a bit of a hack since I would expect the rows to align
    // with out this extra spacing on the first line
    prefix + "[\n  ",
    "\n];"))<<suffix);
}

IGL_INLINE Eigen::IOFormat igl::matlab_format()
{
  // M = [1 2 3;4 5 6];
  // M.format(matlab_format()) produces:
  // [
  //   1 2 3
  //   4 5 6
  // ];
  return Eigen::IOFormat(
    Eigen::FullPrecision,
    0,
    " ",
    "\n",
    "",
    "",
    // This seems like a bit of a hack since I would expect the rows to align
    // with out this extra spacing on the first line
    "[\n  ",
    "\n];");
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciations
template std::basic_string<char, std::char_traits<char>, std::allocator<char> > const igl::matlab_format<double>(Eigen::SparseMatrix<double, 0, int> const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >);
template Eigen::WithFormat<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const igl::matlab_format<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::string);
template Eigen::WithFormat<Eigen::Array<int, -1, -1, 0, -1, -1> > const igl::matlab_format<Eigen::Array<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Array<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >);
template Eigen::WithFormat<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const igl::matlab_format<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >);
#endif
