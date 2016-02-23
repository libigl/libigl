#include "readBF.h"
#include "list_to_matrix.h"
#include <vector>
#include <cstdio>
#include <fstream>
template < 
  typename DerivedWI,
  typename DerivedP,
  typename DerivedC>
IGL_INLINE bool igl::readBF(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedWI> & WI,
  Eigen::PlainObjectBase<DerivedP> & P,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  using namespace std;
  ifstream is(filename);
  if(!is.is_open())
  {
    return false;
  }
  string line;
  std::vector<typename DerivedWI::Scalar> vWI;
  std::vector<typename DerivedP::Scalar> vP;
  std::vector<std::vector<typename DerivedC::Scalar> > vC;
  while(getline(is, line))
  {
    int wi,p;
    double cx,cy,cz;
    if(sscanf(line.c_str(), "%d %d %lg %lg %lg",&wi,&p,&cx,&cy,&cz) != 5)
    {
      return false;
    }
    vWI.push_back(wi);
    vP.push_back(p);
    vC.push_back({cx,cy,cz});
  }
  list_to_matrix(vWI,WI);
  list_to_matrix(vP,P);
  list_to_matrix(vC,C);
  return true;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::readBF<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
