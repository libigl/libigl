#include "igl/matlab/MatlabWorkspace.h"

// IGL
#include "igl/list_to_matrix.h"

// MATLAB
#include "mat.h"

// STL
#include <iostream>
#include <algorithm>

IGL_INLINE igl::MatlabWorkspace::MatlabWorkspace()
{
}

IGL_INLINE igl::MatlabWorkspace::~MatlabWorkspace()
{
  // clean up data
  clear();
}

IGL_INLINE void igl::MatlabWorkspace::clear()
{
  for_each(data.begin(),data.end(),&mxDestroyArray);
}

IGL_INLINE bool igl::MatlabWorkspace::write(const std::string & path) const
{
  using namespace std;
  MATFile * mat_file = matOpen(path.c_str(), "w");
  assert(names.size() == data.size());
  // loop over names and data
  for(int i = 0;i < (int)names.size(); i++)
  {
    // Put variable as LOCAL variable
    int status = matPutVariable(mat_file,names[i].c_str(), data[i]);
    if(status != 0) 
    {
      cerr<<"^MatlabWorkspace::save Error: matPutVariable ("<<names[i]<<
        ") failed"<<endl;
      return false;
    } 
  }
  if(matClose(mat_file) != 0)
  {
    fprintf(stderr,"Error closing file %s\n",path.c_str());
    return false;
  }
  return true;
}

// Treat everything as a double
template <typename DerivedM>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save(
  const Eigen::PlainObjectBase<DerivedM>& M,
  const std::string & name)
{
  using namespace std;
  const int m = M.rows();
  const int n = M.cols();
  mxArray * mx_data = mxCreateDoubleMatrix(m,n,mxREAL);
  data.push_back(mx_data);
  names.push_back(name);
  // Copy data immediately
  // Q: Won't this be incorrect for integers?
  copy(M.data(),M.data()+M.size(),mxGetPr(mx_data));
  return *this;
}

// Treat everything as a double
template <typename MT>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save(
  const Eigen::SparseMatrix<MT>& M,
  const std::string & name)
{
  using namespace std;
  const int m = M.rows();
  const int n = M.cols();
  // THIS WILL NOT WORK FOR ROW-MAJOR
  assert(n==M.outerSize());
  const int nzmax = M.nonZeros();
  mxArray * mx_data = mxCreateSparse(m, n, nzmax, mxREAL);
  data.push_back(mx_data);
  names.push_back(name);
  // Copy data immediately
  double * pr = mxGetPr(mx_data);
  mwIndex * ir = mxGetIr(mx_data);
  mwIndex * jc = mxGetJc(mx_data);

  // Iterate over outside
  int k = 0;
  for(int j=0; j<M.outerSize();j++)
  {
    jc[j] = k;
    // Iterate over inside
    for(typename Eigen::SparseMatrix<MT>::InnerIterator it (M,j); it; ++it)
    {
      pr[k] = it.value();
      ir[k] = it.row();
      k++;
    }
  }
  jc[M.outerSize()] = k;

  return *this;
}

template <typename ScalarM>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save(
  const std::vector<std::vector<ScalarM> > & vM,
  const std::string & name)
{
  Eigen::MatrixXd M;
  list_to_matrix(vM,M);
  return this->save(M,name);
}

template <typename ScalarV>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save(
  const std::vector<ScalarV> & vV,
  const std::string & name)
{
  Eigen::MatrixXd V;
  list_to_matrix(vV,V);
  return this->save(V,name);
}

template <typename DerivedM>
IGL_INLINE igl::MatlabWorkspace& 
  igl::MatlabWorkspace::save_index(
    const Eigen::PlainObjectBase<DerivedM>& M,
    const std::string & name)
{
  DerivedM Mp1 = M;
  Mp1.array() += 1;
  return this->save(Mp1,name);
}

template <typename ScalarM>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save_index(
  const std::vector<std::vector<ScalarM> > & vM,
  const std::string & name)
{
  Eigen::MatrixXd M;
  list_to_matrix(vM,M);
  return this->save_index(M,name);
}

template <typename ScalarV>
IGL_INLINE igl::MatlabWorkspace& igl::MatlabWorkspace::save_index(
  const std::vector<ScalarV> & vV,
  const std::string & name)
{
  Eigen::MatrixXd V;
  list_to_matrix(vV,V);
  return this->save_index(V,name);
}


//template <typename Data>
//bool igl::MatlabWorkspace::save(const Data & M, const std::string & name)
//{
//  using namespace std;
//  // If I don't know the type then I can't save it
//  cerr<<"^MatlabWorkspace::save Error: Unknown data type. "<<
//    name<<" not saved."<<endl;
//  return false;
//}
