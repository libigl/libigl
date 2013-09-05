#ifndef IGL_WRITE_MATLAB_WORKSPACE
#define IGL_WRITE_MATLAB_WORKSPACE
#include "igl/igl_inline.h"

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "mat.h"

namespace igl
{
  // Class which contains data of a matlab workspace which can be written to a
  // .mat file and loaded from matlab
  // 
  // This depends on matlab at compile time (though it shouldn't necessarily
  // have to) but it does not depend on running the matlab engine at run-time.
  //
  // Known bugs: Treats all matrices as doubles (this may actually be desired
  // for some "index" matrices since matlab's sparse command takes doubles
  // rather than int class matrices). It is of course not desired when dealing
  // with logicals or uint's for images.
  class MatlabWorkspace
  {
    private:
      // List of names
      std::vector<std::string> names;
      // List of data pointers
      std::vector<mxArray*> data;
    public:
      MatlabWorkspace();
      ~MatlabWorkspace();
      // Clear names and data of variables in workspace
      IGL_INLINE void clear();
      // Save current list of variables
      //
      // Inputs:
      //   path  path to .mat file
      // Returns true on success, false on failure
      IGL_INLINE bool write(const std::string & path) const;
      // Assign data to a variable name in the workspace
      //
      // Template: 
      //   DerivedM  eigen matrix (e.g. MatrixXd)
      // Inputs:
      //   M  data (usually a matrix)
      //   name  variable name to save into work space
      // Returns true on success, false on failure
      //
      // Known Bugs: Assumes Eigen is using column major ordering
      template <typename DerivedM>
      IGL_INLINE MatlabWorkspace& save(
        const Eigen::PlainObjectBase<DerivedM>& M,
        const std::string & name);
      // Template:
      //   MT  sparse matrix type (e.g. double)
      template <typename MT>
      IGL_INLINE MatlabWorkspace& save(
        const Eigen::SparseMatrix<MT>& M,
        const std::string & name);
      // Templates:
      //   ScalarM  scalar type, e.g. double
      template <typename ScalarM>
      IGL_INLINE MatlabWorkspace& save(
        const std::vector<std::vector<ScalarM> > & vM,
        const std::string & name);
      // Templates:
      //   ScalarV  scalar type, e.g. double
      template <typename ScalarV>
      IGL_INLINE MatlabWorkspace& save(
        const std::vector<ScalarV> & vV,
        const std::string & name);
      // Same as save() but adds 1 to each element, useful for saving "index"
      // matrices like lists of faces or elements
      template <typename DerivedM>
      IGL_INLINE MatlabWorkspace& save_index(
        const Eigen::PlainObjectBase<DerivedM>& M,
        const std::string & name);
      template <typename ScalarM>
      IGL_INLINE MatlabWorkspace& save_index(
        const std::vector<std::vector<ScalarM> > & vM,
        const std::string & name);
      template <typename ScalarV>
      IGL_INLINE MatlabWorkspace& save_index(
        const std::vector<ScalarV> & vV,
        const std::string & name);
  };
}

// Be sure that this is not compiled into libigl.a
#ifdef IGL_HEADER_ONLY
#  include "MatlabWorkspace.cpp"
#endif

#endif

