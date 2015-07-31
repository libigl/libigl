// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/comiso/miq.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>

// includes for VertexIndexing
#include <igl/HalfEdgeIterator.h>
#include <igl/is_border_vertex.h>
#include <igl/vertex_triangle_adjacency.h>


// includes for poissonSolver
#include <gmm/gmm.h>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <CoMISo/Solver/MISolver.hh>
#include <CoMISo/Solver/GMM_Tools.hh>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>

//
#include <igl/cross_field_missmatch.h>
#include <igl/comb_frame_field.h>
#include <igl/comb_cross_field.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/rotate_vectors.h>


// #define DEBUG_PRINT
#include <fstream>
#include <iostream>
#include <igl/matlab_format.h>

#warning "using namespace *; in global scope **must** be removed"
using namespace std;
using namespace Eigen;

#define DEBUGPRINT 0


namespace igl {
namespace comiso {

  class SparseMatrixData{
  protected:
    unsigned int m_nrows;
    unsigned int m_ncols;
    std::vector<unsigned int> m_rowind;
    std::vector<unsigned int> m_colind;
    std::vector<double>       m_vals;

  public:
    unsigned int   nrows()    { return m_nrows      ; }
    unsigned int   ncols()    { return m_ncols      ; }
    unsigned int   nentries() { return m_vals.size(); }
    std::vector<unsigned int>&  rowind()   { return m_rowind     ; }
    std::vector<unsigned int>&  colind()   { return m_colind     ; }
    std::vector<double>&        vals()     { return m_vals       ; }

    // create an empty matrix with a fixed number of rows
    IGL_INLINE SparseMatrixData()
    {
      initialize(0,0);
    }

    // create an empty matrix with a fixed number of rows
    IGL_INLINE void initialize(int nr, int nc) {
      assert(nr >= 0 && nc >=0);
      m_nrows = nr;
      m_ncols = nc;

      m_rowind.resize(0);
      m_colind.resize(0);
      m_vals.resize(0);
    }

    // add a nonzero entry to the matrix
    // no checks are done for coinciding entries
    // the interpretation of the repeated entries (replace or add)
    // depends on how the actual sparse matrix datastructure is constructed

    IGL_INLINE void addEntryCmplx(unsigned int i, unsigned int j, std::complex<double> val) {
      m_rowind.push_back(2*i);   m_colind.push_back(2*j);   m_vals.push_back( val.real());
      m_rowind.push_back(2*i);   m_colind.push_back(2*j+1); m_vals.push_back(-val.imag());
      m_rowind.push_back(2*i+1); m_colind.push_back(2*j);   m_vals.push_back( val.imag());
      m_rowind.push_back(2*i+1); m_colind.push_back(2*j+1); m_vals.push_back( val.real());
    }

    IGL_INLINE void addEntryReal(unsigned int i, unsigned int j, double val) {
      m_rowind.push_back(i);   m_colind.push_back(j);   m_vals.push_back(val);
    }

    IGL_INLINE virtual ~SparseMatrixData() {
    }

  };

  // a small class to manage storage for matrix data
  // not using stl vectors: want to make all memory management
  // explicit to avoid hidden automatic reallocation
  // TODO: redo with STL vectors but with explicit mem. management

  class SparseSystemData {
  private:
    // matrix representation,  A[rowind[i],colind[i]] = vals[i]
    // right-hand side
    SparseMatrixData m_A;
    double       *m_b;
    double       *m_x;

  public:
    IGL_INLINE SparseMatrixData& A() { return m_A; }
    IGL_INLINE double*        b()        { return m_b       ; }
    IGL_INLINE double*        x()        { return m_x       ; }
    IGL_INLINE unsigned int   nrows()    { return  m_A.nrows(); }

  public:

    IGL_INLINE SparseSystemData(): m_A(), m_b(NULL), m_x(NULL){ }

    IGL_INLINE void initialize(unsigned int nr, unsigned int nc) {
      m_A.initialize(nr,nc);
      m_b      = new          double[nr];
      m_x      = new          double[nr];
      assert(m_b);
      std::fill( m_b,  m_b+nr, 0.);
    }

    IGL_INLINE void addRHSCmplx(unsigned int i, std::complex<double> val) {
      assert( 2*i+1 < m_A.nrows());
      m_b[2*i] += val.real(); m_b[2*i+1] += val.imag();
    }

    IGL_INLINE void setRHSCmplx(unsigned int i, std::complex<double> val) {
      assert( 2*i+1 < m_A.nrows());
      m_b[2*i] = val.real(); m_b[2*i+1] = val.imag();
    }

    IGL_INLINE std::complex<double> getRHSCmplx(unsigned int i) {
      assert( 2*i+1 < m_A.nrows());
      return std::complex<double>( m_b[2*i], m_b[2*i+1]);
    }

    IGL_INLINE double getRHSReal(unsigned int i) {
      assert( i < m_A.nrows());
      return m_b[i];
    }

    IGL_INLINE std::complex<double> getXCmplx(unsigned int i) {
      assert( 2*i+1 < m_A.nrows());
      return std::complex<double>( m_x[2*i], m_x[2*i+1]);
    }

    IGL_INLINE void cleanMem() {
      //m_A.cleanup();
      delete [] m_b;
      delete [] m_x;
    }

    IGL_INLINE virtual ~SparseSystemData() {
      delete [] m_b;
      delete [] m_x;
    }
  };

  struct SeamInfo
  {
    int v0,v0p,v1,v1p;
    int integerVar;
    unsigned char MMatch;

    IGL_INLINE SeamInfo(int _v0,
                        int _v1,
                        int _v0p,
                        int _v1p,
                        int _MMatch,
                        int _integerVar);

    IGL_INLINE SeamInfo(const SeamInfo &S1);
  };

  struct MeshSystemInfo
  {
    ///total number of scalar variables
    int num_scalar_variables;
    ////number of vertices variables
    int num_vert_variables;
    ///num of integer for cuts
    int num_integer_cuts;
    ///this are used for drawing purposes
    std::vector<SeamInfo> EdgeSeamInfo;
#if 0
    ///this are values of integer variables after optimization
    std::vector<int> IntegerValues;
#endif
  };


  template <typename DerivedV, typename DerivedF>
  class VertexIndexing
  {
  public:
    // Input:
    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F;
    const Eigen::PlainObjectBase<DerivedF> &TT;
    const Eigen::PlainObjectBase<DerivedF> &TTi;
    // const Eigen::PlainObjectBase<DerivedV> &PD1;
    // const Eigen::PlainObjectBase<DerivedV> &PD2;

    const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_MMatch;
    // const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular; // bool
    // const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_SingularDegree; // vertex;
    const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_Seams; // 3 bool


    ///this handle for mesh TODO: move with the other global variables
    MeshSystemInfo Handle_SystemInfo;

    // Output:
    ///this maps the integer for edge - face
    Eigen::MatrixXi Handle_Integer; // TODO: remove it is useless

    ///per face indexes of vertex in the solver
    Eigen::MatrixXi HandleS_Index;

    ///per vertex variable indexes
    std::vector<std::vector<int> > HandleV_Integer;

    // internal
    std::vector<std::vector<int> > VF, VFi;
    std::vector<bool> V_border; // bool

    IGL_INLINE VertexIndexing(const Eigen::PlainObjectBase<DerivedV> &_V,
                              const Eigen::PlainObjectBase<DerivedF> &_F,
                              const Eigen::PlainObjectBase<DerivedF> &_TT,
                              const Eigen::PlainObjectBase<DerivedF> &_TTi,
                              //  const Eigen::PlainObjectBase<DerivedV> &_PD1,
                              //  const Eigen::PlainObjectBase<DerivedV> &_PD2,
                              const Eigen::Matrix<int, Eigen::Dynamic, 3> &_Handle_MMatch,
                              //  const Eigen::Matrix<int, Eigen::Dynamic, 1> &_Handle_Singular,
                              //  const Eigen::Matrix<int, Eigen::Dynamic, 1> &_Handle_SingularDegree,
                              const Eigen::Matrix<int, Eigen::Dynamic, 3> &_Handle_Seams
                              );

    ///vertex to variable mapping
    IGL_INLINE void InitMapping();

    IGL_INLINE void InitFaceIntegerVal();

    IGL_INLINE void InitSeamInfo();


  private:
    ///this maps back index to vertices
    std::vector<int> IndexToVert; // TODO remove it is useless

    ///this is used for drawing purposes
    std::vector<int> duplicated; // TODO remove it is useless

    IGL_INLINE void FirstPos(const int v, int &f, int &edge);

    IGL_INLINE int AddNewIndex(const int v0);

    IGL_INLINE bool HasIndex(int indexVert,int indexVar);

    IGL_INLINE void GetSeamInfo(const int f0,
                                const int f1,
                                const int indexE,
                                int &v0,int &v1,
                                int &v0p,int &v1p,
                                unsigned char &_MMatch,
                                int &integerVar);
    IGL_INLINE bool IsSeam(const int f0, const int f1);

    ///find initial position of the pos to
    // assing face to vert inxex correctly
    IGL_INLINE void FindInitialPos(const int vert, int &edge, int &face);


    ///intialize the mapping given an initial pos
    ///whih must be initialized with FindInitialPos
    IGL_INLINE void MapIndexes(const int  vert, const int edge_init, const int f_init);

    ///intialize the mapping for a given vertex
    IGL_INLINE void InitMappingSeam(const int vert);

    ///intialize the mapping for a given sampled mesh
    IGL_INLINE void InitMappingSeam();

    ///test consistency of face variables per vert mapping
    IGL_INLINE void TestSeamMappingFace(const int f);

    ///test consistency of face variables per vert mapping
    IGL_INLINE void TestSeamMappingVertex(int indexVert);

    ///check consistency of variable mapping across seams
    IGL_INLINE void TestSeamMapping();

  };


  template <typename DerivedV, typename DerivedF>
  class PoissonSolver
  {

  public:
    IGL_INLINE void SolvePoisson(Eigen::VectorXd Stiffness,
                                 double vector_field_scale=0.1f,
                                 double grid_res=1.f,
                                 bool direct_round=true,
                                 int localIter=0,
                                 bool _integer_rounding=true,
                                 bool _singularity_rounding=true,
                                 std::vector<int> roundVertices = std::vector<int>(),
                                 std::vector<std::vector<int> > hardFeatures = std::vector<std::vector<int> >());

    IGL_INLINE PoissonSolver(const Eigen::PlainObjectBase<DerivedV> &_V,
                             const Eigen::PlainObjectBase<DerivedF> &_F,
                             const Eigen::PlainObjectBase<DerivedF> &_TT,
                             const Eigen::PlainObjectBase<DerivedF> &_TTi,
                             const Eigen::PlainObjectBase<DerivedV> &_PD1,
                             const Eigen::PlainObjectBase<DerivedV> &_PD2,
                             const Eigen::MatrixXi &_HandleS_Index,
                             const Eigen::Matrix<int, Eigen::Dynamic, 1>&_Handle_Singular,
                             const MeshSystemInfo &_Handle_SystemInfo
                             );

    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F;
    const Eigen::PlainObjectBase<DerivedF> &TT;
    const Eigen::PlainObjectBase<DerivedF> &TTi;
    const Eigen::PlainObjectBase<DerivedV> &PD1;
    const Eigen::PlainObjectBase<DerivedV> &PD2;
    const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular; // bool
    const Eigen::MatrixXi &HandleS_Index; //todo

    const MeshSystemInfo &Handle_SystemInfo;

    // Internal:
    Eigen::MatrixXd doublearea;
    Eigen::VectorXd Handle_Stiffness;
    Eigen::PlainObjectBase<DerivedV> N;
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFi;
    Eigen::MatrixXd UV; // this is probably useless

    // Output:
    // per wedge UV coordinates, 6 coordinates (1 face) per row
    Eigen::MatrixXd WUV;

    ///solver data
    SparseSystemData S;

    ///vector of unknowns
    std::vector< double > X;

    ////REAL PART
    ///number of fixed vertex
    unsigned int n_fixed_vars;

    ///the number of REAL variables for vertices
    unsigned int n_vert_vars;

    ///total number of variables of the system,
    ///do not consider constraints, but consider integer vars
    unsigned int num_total_vars;

    //////INTEGER PART
    ///the total number of integer variables
    unsigned int n_integer_vars;

    ///CONSTRAINT PART
    ///number of cuts constraints
    unsigned int num_cut_constraint;

    // number of user-defined constraints
    unsigned int num_userdefined_constraint;

    ///total number of constraints equations
    unsigned int num_constraint_equations;

    ///total size of the system including constraints
    unsigned int system_size;

    ///if you intend to make integer rotation
    ///and translations
    bool integer_jumps_bary;

    ///vector of blocked vertices
    std::vector<int> Hard_constraints;

    ///vector of indexes to round
    std::vector<int> ids_to_round;

    ///vector of indexes to round
    std::vector<std::vector<int > > userdefined_constraints;

    ///boolean that is true if rounding to integer is needed
    bool integer_rounding;

    ///START SYSTEM ACCESS METHODS
    ///add an entry to the LHS
    IGL_INLINE void AddValA(int Xindex,
                            int Yindex,
                            double val);

    ///add a complex entry to the LHS
    IGL_INLINE void AddComplexA(int VarXindex,
                                int VarYindex,
                                std::complex<double> val);

    ///add a velue to the RHS
    IGL_INLINE void AddValB(int Xindex,
                            double val);

    ///add the area term, scalefactor is used to sum up
    ///and normalize on the overlap zones
    IGL_INLINE void AddAreaTerm(int index[3][3][2],double ScaleFactor);

    ///set the diagonal of the matrix (which is zero at the beginning)
    ///such that the sum of a row or a colums is zero
    IGL_INLINE void SetDiagonal(double val[3][3]);

    ///given a vector of scalar values and
    ///a vector of indexes add such values
    ///as specified by the indexes
    IGL_INLINE void AddRHS(double b[6],
                           int index[3]);

    ///add a 3x3 block matrix to the system matrix...
    ///indexes are specified in the 3x3 matrix of x,y pairs
    ///indexes must be multiplied by 2 cause u and v
    IGL_INLINE void Add33Block(double val[3][3], int index[3][3][2]);

    ///add a 3x3 block matrix to the system matrix...
    ///indexes are specified in the 3x3 matrix of x,y pairs
    ///indexes must be multiplied by 2 cause u and v
    IGL_INLINE void Add44Block(double val[4][4],int index[4][4][2]);
    ///END SYSTEM ACCESS METHODS

    ///START COMMON MATH FUNCTIONS
    ///return the complex encoding the rotation
    ///for a given missmatch interval
    IGL_INLINE std::complex<double> GetRotationComplex(int interval);
    ///END COMMON MATH FUNCTIONS


    ///START ENERGY MINIMIZATION PART
    ///initialize the LHS for a given face
    ///for minimization of Dirichlet's energy
    IGL_INLINE void perElementLHS(int f,
                                  double val[3][3],
                                  int index[3][3][2]);

    ///initialize the RHS for a given face
    ///for minimization of Dirichlet's energy
    IGL_INLINE void perElementRHS(int f,
                                  double b[6],
                                  double vector_field_scale=1);

    ///evaluate the LHS and RHS for a single face
    ///for minimization of Dirichlet's energy
    IGL_INLINE void PerElementSystemReal(int f,
                                         double val[3][3],
                                         int index[3][3][2],
                                         double b[6],
                                         double vector_field_scale=1.0);
    ///END ENERGY MINIMIZATION PART

    ///START FIXING VERTICES
    ///set a given vertex as fixed
    IGL_INLINE void AddFixedVertex(int v);

    ///find vertex to fix in case we're using
    ///a vector field NB: multiple components not handled
    IGL_INLINE void FindFixedVertField();

    ///find hard constraint depending if using or not
    ///a vector field
    IGL_INLINE void FindFixedVert();

    IGL_INLINE int GetFirstVertexIndex(int v);

    ///fix the vertices which are flagged as fixed
    IGL_INLINE void FixBlockedVertex();
    ///END FIXING VERTICES

    ///HANDLING SINGULARITY
    //set the singularity round to integer location
    IGL_INLINE void AddSingularityRound();

    IGL_INLINE void AddToRoundVertices(std::vector<int> ids);

    ///START GENERIC SYSTEM FUNCTIONS
    //build the laplacian matrix cyclyng over all rangemaps
    //and over all faces
    IGL_INLINE void BuildLaplacianMatrix(double vfscale=1);

    ///find different sized of the system
    IGL_INLINE void FindSizes();

    IGL_INLINE void AllocateSystem();

    ///intitialize the whole matrix
    IGL_INLINE void InitMatrix();

    ///map back coordinates after that
    ///the system has been solved
    IGL_INLINE void MapCoords();
    ///END GENERIC SYSTEM FUNCTIONS

    ///set the constraints for the inter-range cuts
    IGL_INLINE void BuildSeamConstraintsExplicitTranslation();

    ///set the constraints for the inter-range cuts
    IGL_INLINE void BuildUserDefinedConstraints();

    ///call of the mixed integer solver
    IGL_INLINE void MixedIntegerSolve(double cone_grid_res=1,
                                      bool direct_round=true,
                                      int localIter=0);

    IGL_INLINE void clearUserConstraint();

    IGL_INLINE void addSharpEdgeConstraint(int fid, int vid);

  };

  template <typename DerivedV, typename DerivedF, typename DerivedU>
  class MIQ_class
  {
  private:
    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F;
    Eigen::MatrixXd WUV;
    // internal
    Eigen::PlainObjectBase<DerivedF> TT;
    Eigen::PlainObjectBase<DerivedF> TTi;

    // Stiffness per face
    Eigen::VectorXd Handle_Stiffness;
    Eigen::PlainObjectBase<DerivedV> B1, B2, B3;

  public:
    IGL_INLINE MIQ_class(const Eigen::PlainObjectBase<DerivedV> &V_,
                         const Eigen::PlainObjectBase<DerivedF> &F_,
                         const Eigen::PlainObjectBase<DerivedV> &PD1_combed,
                         const Eigen::PlainObjectBase<DerivedV> &PD2_combed,
                         // const Eigen::PlainObjectBase<DerivedV> &BIS1_combed,
                         // const Eigen::PlainObjectBase<DerivedV> &BIS2_combed,
                         const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_MMatch,
                         const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular,
                         // const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_SingularDegree,
                         const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_Seams,
                         Eigen::PlainObjectBase<DerivedU> &UV,
                         Eigen::PlainObjectBase<DerivedF> &FUV,
                         double GradientSize = 30.0,
                         double Stiffness = 5.0,
                         bool DirectRound = false,
                         int iter = 5,
                         int localIter = 5,
                         bool DoRound = true,
                         bool SingularityRound=true,
                         std::vector<int> roundVertices = std::vector<int>(),
                         std::vector<std::vector<int> > hardFeatures = std::vector<std::vector<int> >());


    IGL_INLINE void extractUV(Eigen::PlainObjectBase<DerivedU> &UV_out,
                              Eigen::PlainObjectBase<DerivedF> &FUV_out);

  private:
    IGL_INLINE int NumFlips(const Eigen::MatrixXd& WUV);

    IGL_INLINE double Distortion(int f, double h, const Eigen::MatrixXd& WUV);

    IGL_INLINE double LaplaceDistortion(const int f, double h, const Eigen::MatrixXd& WUV);

    IGL_INLINE bool updateStiffeningJacobianDistorsion(double grad_size, const Eigen::MatrixXd& WUV);

    IGL_INLINE bool IsFlipped(const Eigen::Vector2d &uv0,
                              const Eigen::Vector2d &uv1,
                              const Eigen::Vector2d &uv2);

    IGL_INLINE bool IsFlipped(const int i, const Eigen::MatrixXd& WUV);

  };
};
}

IGL_INLINE igl::comiso::SeamInfo::SeamInfo(int _v0,
                                   int _v1,
                                   int _v0p,
                                   int _v1p,
                                   int _MMatch,
                                   int _integerVar)
{
  v0=_v0;
  v1=_v1;
  v0p=_v0p;
  v1p=_v1p;
  integerVar=_integerVar;
  MMatch=_MMatch;
}

IGL_INLINE igl::comiso::SeamInfo::SeamInfo(const SeamInfo &S1)
{
  v0=S1.v0;
  v1=S1.v1;
  v0p=S1.v0p;
  v1p=S1.v1p;
  integerVar=S1.integerVar;
  MMatch=S1.MMatch;
}


template <typename DerivedV, typename DerivedF>
IGL_INLINE igl::comiso::VertexIndexing<DerivedV, DerivedF>::VertexIndexing(const Eigen::PlainObjectBase<DerivedV> &_V,
                                                                   const Eigen::PlainObjectBase<DerivedF> &_F,
                                                                   const Eigen::PlainObjectBase<DerivedF> &_TT,
                                                                   const Eigen::PlainObjectBase<DerivedF> &_TTi,
                                                                   // const Eigen::PlainObjectBase<DerivedV> &_PD1,
                                                                   // const Eigen::PlainObjectBase<DerivedV> &_PD2,
                                                                   const Eigen::Matrix<int, Eigen::Dynamic, 3> &_Handle_MMatch,
                                                                   // const Eigen::Matrix<int, Eigen::Dynamic, 1> &_Handle_Singular,
                                                                   // const Eigen::Matrix<int, Eigen::Dynamic, 1> &_Handle_SingularDegree,
                                                                   const Eigen::Matrix<int, Eigen::Dynamic, 3> &_Handle_Seams

                                                                   ):
V(_V),
F(_F),
TT(_TT),
TTi(_TTi),
// PD1(_PD1),
// PD2(_PD2),
Handle_MMatch(_Handle_MMatch),
// Handle_Singular(_Handle_Singular),
// Handle_SingularDegree(_Handle_SingularDegree),
Handle_Seams(_Handle_Seams)
{
  #ifdef DEBUG_PRINT
  cerr<<igl::matlab_format(Handle_Seams,"Handle_Seams");
#endif
  V_border = igl::is_border_vertex(V,F);
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  IndexToVert.clear();

  Handle_SystemInfo.num_scalar_variables=0;
  Handle_SystemInfo.num_vert_variables=0;
  Handle_SystemInfo.num_integer_cuts=0;

  duplicated.clear();

  HandleS_Index = Eigen::MatrixXi::Constant(F.rows(),3,-1);

  Handle_Integer = Eigen::MatrixXi::Constant(F.rows(),3,-1);

  HandleV_Integer.resize(V.rows());
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::FirstPos(const int v, int &f, int &edge)
{
  f    = VF[v][0];  // f=v->cVFp();
  edge = VFi[v][0]; // edge=v->cVFi();
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE int igl::comiso::VertexIndexing<DerivedV, DerivedF>::AddNewIndex(const int v0)
{
  Handle_SystemInfo.num_scalar_variables++;
  HandleV_Integer[v0].push_back(Handle_SystemInfo.num_scalar_variables);
  IndexToVert.push_back(v0);
  return Handle_SystemInfo.num_scalar_variables;
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::comiso::VertexIndexing<DerivedV, DerivedF>::HasIndex(int indexVert,int indexVar)
{
  for (unsigned int i=0;i<HandleV_Integer[indexVert].size();i++)
    if (HandleV_Integer[indexVert][i]==indexVar)return true;
  return false;
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::GetSeamInfo(const int f0,
                                                                     const int f1,
                                                                     const int indexE,
                                                                     int &v0,int &v1,
                                                                     int &v0p,int &v1p,
                                                                     unsigned char &_MMatch,
                                                                     int &integerVar)
{
  int edgef0 = indexE;
  v0 = HandleS_Index(f0,edgef0);
  v1 = HandleS_Index(f0,(edgef0+1)%3);
  ////get the index on opposite side
  assert(TT(f0,edgef0) == f1);
  int edgef1 = TTi(f0,edgef0);
  v1p = HandleS_Index(f1,edgef1);
  v0p = HandleS_Index(f1,(edgef1+1)%3);

  integerVar = Handle_Integer(f0,edgef0);
  _MMatch = Handle_MMatch(f0,edgef0);
  assert(F(f0,edgef0)         == F(f1,((edgef1+1)%3)));
  assert(F(f0,((edgef0+1)%3)) == F(f1,edgef1));
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::comiso::VertexIndexing<DerivedV, DerivedF>::IsSeam(const int f0, const int f1)
{
  for (int i=0;i<3;i++)
  {
    int f_clos = TT(f0,i);

    if (f_clos == -1)
      continue; ///border

    if (f_clos == f1)
      return(Handle_Seams(f0,i));
  }
  assert(0);
  return false;
}

///find initial position of the pos to
// assing face to vert inxex correctly
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::FindInitialPos(const int vert,
                                                                        int &edge,
                                                                        int &face)
{
  int f_init;
  int edge_init;
  FirstPos(vert,f_init,edge_init); // todo manually IGL_INLINE the function
  igl::HalfEdgeIterator<DerivedF> VFI(F,TT,TTi,f_init,edge_init);

#ifdef DEBUG_PRINT
  cerr<<"--FindInitialPos--"<<endl;
#endif
  bool vertexB = V_border[vert];
  bool possible_split=false;
  bool complete_turn=false;
  do
  {
    int curr_f = VFI.Fi();
    int curr_edge=VFI.Ei();
#ifdef DEBUG_PRINT
    cerr<<"@ face "<<curr_f<<", edge "<< F(curr_f,curr_edge)<<" - "<< F(curr_f,(curr_edge+1)%3)<<endl;
#endif
    VFI.NextFE();
    int next_f=VFI.Fi();
#ifdef DEBUG_PRINT
    cerr<<"next face "<<next_f<<", edge "<< F(next_f,VFI.Ei())<<" - "<< F(next_f,(VFI.Ei()+1)%3)<<endl;
#endif
    ///test if I've just crossed a border
    bool on_border=(TT(curr_f,curr_edge)==-1);
#ifdef DEBUG_PRINT
    cerr<<"on_border: "<<on_border<<endl;
#endif
    //bool mismatch=false;
    bool seam=false;

    #ifdef DEBUG_PRINT
    cerr<<igl::matlab_format(Handle_Seams,"Handle_Seams");
    #endif
    ///or if I've just crossed a seam
    ///if I'm on a border I MUST start from the one next t othe border
    if (!vertexB)
      //seam=curr_f->IsSeam(next_f);
      seam=IsSeam(curr_f,next_f);
    // if (vertexB)
    // assert(!Handle_Singular(vert));
    // ;
    //assert(!vert->IsSingular());
#ifdef DEBUG_PRINT
    cerr<<"seam: "<<seam<<endl;
#endif
    possible_split=((on_border)||(seam));
#ifdef DEBUG_PRINT
    cerr<<"possible_split: "<<possible_split<<endl;
#endif
    complete_turn = next_f == f_init;
#ifdef DEBUG_PRINT
    cerr<<"complete_turn: "<<complete_turn<<endl;
#endif
  } while ((!possible_split)&&(!complete_turn));
  face=VFI.Fi();
  edge=VFI.Ei();
#ifdef DEBUG_PRINT
  cerr<<"FindInitialPos done. Face: "<<face<<", edge: "<< F(face,edge)<<" - "<< F(face,(edge+1)%3)<<endl;
#endif
  ///test that is not on a border
  //assert(face->FFp(edge)!=face);
}



///intialize the mapping given an initial pos
///whih must be initialized with FindInitialPos
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::MapIndexes(const int  vert,
                                                                    const int edge_init,
                                                                    const int f_init)
{
  ///check that is not on border..
  ///in such case maybe it's non manyfold
  ///insert an initial index
  int curr_index=AddNewIndex(vert);
#ifdef DEBUG_PRINT
  cerr<<"--MapIndexes--"<<endl;
#endif
#ifdef DEBUG_PRINT
  cerr<<"adding vertex for "<<vert<<endl;
#endif
  ///and initialize the jumping pos
  igl::HalfEdgeIterator<DerivedF> VFI(F,TT,TTi,f_init,edge_init);
  bool complete_turn=false;
  do
  {
    int curr_f = VFI.Fi();
    int curr_edge = VFI.Ei();
#ifdef DEBUG_PRINT
    cerr<<"Adding vertex "<<curr_index<<" to face "<<curr_f<<", edge "<< F(curr_f,curr_edge)<<" - "<< F(curr_f,(curr_edge+1)%3)<<endl;
#endif
    ///assing the current index
    HandleS_Index(curr_f,curr_edge) = curr_index;
#ifdef DEBUG_PRINT
    cerr<<igl::matlab_format(HandleS_Index,"HandleS_Index")<<endl;
#endif
    VFI.NextFE();
    int next_f = VFI.Fi();
#ifdef DEBUG_PRINT
    cerr<<"next face "<<next_f<<", edge "<< F(next_f,VFI.Ei())<<" - "<< F(next_f,(VFI.Ei()+1)%3)<<endl;
#endif
    ///test if I've finiseh with the face exploration
    complete_turn = (next_f==f_init);
#ifdef DEBUG_PRINT
    cerr<<"complete_turn: "<<complete_turn<<endl;
#endif
    ///or if I've just crossed a mismatch
    if (!complete_turn)
    {
      bool seam=false;
      //seam=curr_f->IsSeam(next_f);
      seam=IsSeam(curr_f,next_f);
      if (seam)
      {
        ///then add a new index
        curr_index=AddNewIndex(vert);
#ifdef DEBUG_PRINT
        cerr<<"Found a seam, adding vertex for "<<vert<<endl;
#endif
      }
    }
  } while (!complete_turn);
}

///intialize the mapping for a given vertex
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::InitMappingSeam(const int vert)
{
  ///first rotate until find the first pos after a mismatch
  ///or a border or return to the first position...
  int f_init = VF[vert][0];
  int indexE = VFi[vert][0];

  igl::HalfEdgeIterator<DerivedF> VFI(F,TT,TTi,f_init,indexE);

  int edge_init;
  int face_init;
#ifdef DEBUG_PRINT
  cerr<<"---Vertex: "<<vert<<"---"<<endl;
#endif
  FindInitialPos(vert,edge_init,face_init);
  MapIndexes(vert,edge_init,face_init);
}

///intialize the mapping for a given sampled mesh
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::InitMappingSeam()
{
  //num_scalar_variables=-1;
  Handle_SystemInfo.num_scalar_variables=-1;
  for (unsigned int i=0;i<V.rows();i++)
    InitMappingSeam(i);

  for (unsigned int j=0;j<V.rows();j++)
  {
    assert(HandleV_Integer[j].size()>0);
    if (HandleV_Integer[j].size()>1)
      duplicated.push_back(j);
  }
}

///test consistency of face variables per vert mapping
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::TestSeamMappingFace(const int f)
{
  for (int k=0;k<3;k++)
  {
    int indexV=HandleS_Index(f,k);
    int v = F(f,k);
    bool has_index=HasIndex(v,indexV);
    assert(has_index);
  }
}

///test consistency of face variables per vert mapping
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::TestSeamMappingVertex(int indexVert)
{
  for (unsigned int k=0;k<HandleV_Integer[indexVert].size();k++)
  {
    int indexV=HandleV_Integer[indexVert][k];

    ///get faces sharing vertex
    std::vector<int> faces = VF[indexVert];
    std::vector<int> indexes = VFi[indexVert];

    for (unsigned int j=0;j<faces.size();j++)
    {
      int f = faces[j];
      int index = indexes[j];
      assert(F(f,index) == indexVert);
      assert((index>=0)&&(index<3));

      if (HandleS_Index(f,index) == indexV)
        return;
    }
  }
  assert(0);
}


///check consistency of variable mapping across seams
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::TestSeamMapping()
{
  printf("\n TESTING SEAM INDEXES \n");
  ///test F-V mapping
  for (unsigned int j=0;j<F.rows();j++)
    TestSeamMappingFace(j);

  ///TEST  V-F MAPPING
  for (unsigned int j=0;j<V.rows();j++)
    TestSeamMappingVertex(j);

}


///vertex to variable mapping
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::InitMapping()
{
  //use_direction_field=_use_direction_field;

  IndexToVert.clear();
  duplicated.clear();

  InitMappingSeam();

  Handle_SystemInfo.num_vert_variables=Handle_SystemInfo.num_scalar_variables+1;

  ///end testing...
  TestSeamMapping();
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::InitFaceIntegerVal()
{
  Handle_SystemInfo.num_integer_cuts=0;
  for (unsigned int j=0;j<F.rows();j++)
  {
    for (int k=0;k<3;k++)
    {
      if (Handle_Seams(j,k))
      {
        Handle_Integer(j,k) = Handle_SystemInfo.num_integer_cuts;
        Handle_SystemInfo.num_integer_cuts++;
      }
      else
        Handle_Integer(j,k)=-1;
    }
  }
}


template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::VertexIndexing<DerivedV, DerivedF>::InitSeamInfo()
{
  Handle_SystemInfo.EdgeSeamInfo.clear();
  for (unsigned int f0=0;f0<F.rows();f0++)
  {
    for (int k=0;k<3;k++)
    {
      int f1 = TT(f0,k);

      if (f1 == -1)
        continue;

      bool seam = Handle_Seams(f0,k);
      if (seam)
      {
        int v0,v0p,v1,v1p;
        unsigned char MM;
        int integerVar;
        GetSeamInfo(f0,f1,k,v0,v1,v0p,v1p,MM,integerVar);
        Handle_SystemInfo.EdgeSeamInfo.push_back(SeamInfo(v0,v1,v0p,v1p,MM,integerVar));
      }
    }
  }
}



template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::SolvePoisson(Eigen::VectorXd Stiffness,
                                                                     double vector_field_scale,
                                                                     double grid_res,
                                                                     bool direct_round,
                                                                     int localIter,
                                                                     bool _integer_rounding,
                                                                     bool _singularity_rounding,
                                                                     std::vector<int> roundVertices,
                                                                     std::vector<std::vector<int> > hardFeatures)
{
  Handle_Stiffness = Stiffness;

  //initialization of flags and data structures
  integer_rounding=_integer_rounding;

  ids_to_round.clear();

  clearUserConstraint();
  // copy the user constraints number
  for (size_t i = 0; i < hardFeatures.size(); ++i)
  {
    addSharpEdgeConstraint(hardFeatures[i][0],hardFeatures[i][1]);
  }

  ///Initializing Matrix

  int t0=clock();

  ///initialize the matrix ALLOCATING SPACE
  InitMatrix();
  if (DEBUGPRINT)
    printf("\n ALLOCATED THE MATRIX \n");

  ///build the laplacian system
  BuildLaplacianMatrix(vector_field_scale);

  // add seam constraints
  BuildSeamConstraintsExplicitTranslation();

  // add user defined constraints
  BuildUserDefinedConstraints();

  ////add the lagrange multiplier
  FixBlockedVertex();

  if (DEBUGPRINT)
    printf("\n BUILT THE MATRIX \n");

  if (integer_rounding)
    AddToRoundVertices(roundVertices);

  if (_singularity_rounding)
    AddSingularityRound();

  int t1=clock();
  if (DEBUGPRINT) printf("\n time:%d \n",t1-t0);
  if (DEBUGPRINT) printf("\n SOLVING \n");

  MixedIntegerSolve(grid_res,direct_round,localIter);

  int t2=clock();
  if (DEBUGPRINT) printf("\n time:%d \n",t2-t1);
  if (DEBUGPRINT) printf("\n ASSIGNING COORDS \n");

  MapCoords();

  int t3=clock();
  if (DEBUGPRINT) printf("\n time:%d \n",t3-t2);
  if (DEBUGPRINT) printf("\n FINISHED \n");
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE igl::comiso::PoissonSolver<DerivedV, DerivedF>
::PoissonSolver(const Eigen::PlainObjectBase<DerivedV> &_V,
                const Eigen::PlainObjectBase<DerivedF> &_F,
                const Eigen::PlainObjectBase<DerivedF> &_TT,
                const Eigen::PlainObjectBase<DerivedF> &_TTi,
                const Eigen::PlainObjectBase<DerivedV> &_PD1,
                const Eigen::PlainObjectBase<DerivedV> &_PD2,
                const Eigen::MatrixXi &_HandleS_Index,
                const Eigen::Matrix<int, Eigen::Dynamic, 1>&_Handle_Singular,
                const MeshSystemInfo &_Handle_SystemInfo //todo: const?
):
V(_V),
F(_F),
TT(_TT),
TTi(_TTi),
PD1(_PD1),
PD2(_PD2),
HandleS_Index(_HandleS_Index),
Handle_Singular(_Handle_Singular),
Handle_SystemInfo(_Handle_SystemInfo)
{
  UV        = Eigen::MatrixXd(V.rows(),2);
  WUV       = Eigen::MatrixXd(F.rows(),6);
  igl::doublearea(V,F,doublearea);
  igl::per_face_normals(V,F,N);
  igl::vertex_triangle_adjacency(V,F,VF,VFi);
}


///START SYSTEM ACCESS METHODS
///add an entry to the LHS
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddValA(int Xindex,
                                                                int Yindex,
                                                                double val)
{
  int size=(int)S.nrows();
  assert(0 <= Xindex && Xindex < size);
  assert(0 <= Yindex && Yindex < size);
  S.A().addEntryReal(Xindex,Yindex,val);
}

///add a complex entry to the LHS
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddComplexA(int VarXindex,
                                                                    int VarYindex,
                                                                    std::complex<double> val)
{
  int size=(int)S.nrows()/2;
  assert(0 <= VarXindex && VarXindex < size);
  assert(0 <= VarYindex && VarYindex < size);
  S.A().addEntryCmplx(VarXindex,VarYindex,val);
}

///add a velue to the RHS
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddValB(int Xindex,
                                                                double val)
{
  int size=(int)S.nrows();
  assert(0 <= Xindex && Xindex < size);
  S.b()[Xindex] += val;
}

///add the area term, scalefactor is used to sum up
///and normalize on the overlap zones
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddAreaTerm(int index[3][3][2],double ScaleFactor)
{
  const double entry = 0.5*ScaleFactor;
  double val[3][3]= {
    {0,       entry, -entry},
    {-entry,      0,  entry},
    {entry,  -entry,      0}
  };

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    {
      ///add for both u and v
      int Xindex=index[i][j][0]*2;
      int Yindex=index[i][j][1]*2;

      AddValA(Xindex+1,Yindex,-val[i][j]);
      AddValA(Xindex,Yindex+1,val[i][j]);
    }
}

///set the diagonal of the matrix (which is zero at the beginning)
///such that the sum of a row or a colums is zero
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::SetDiagonal(double val[3][3])
{
  for (int i=0;i<3;i++)
  {
    double sum=0;
    for (int j=0;j<3;j++)
      sum+=val[i][j];
    val[i][i]=-sum;
  }
}

///given a vector of scalar values and
///a vector of indexes add such values
///as specified by the indexes
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddRHS(double b[6],
                                                               int index[3])
{
  for (int i=0;i<3;i++)
  {
    double valU=b[i*2];
    double valV=b[(i*2)+1];
    AddValB((index[i]*2),valU);
    AddValB((index[i]*2)+1,valV);
  }
}

///add a 3x3 block matrix to the system matrix...
///indexes are specified in the 3x3 matrix of x,y pairs
///indexes must be multiplied by 2 cause u and v
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::Add33Block(double val[3][3], int index[3][3][2])
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    {
      ///add for both u and v
      int Xindex=index[i][j][0]*2;
      int Yindex=index[i][j][1]*2;
      assert((unsigned)Xindex<(n_vert_vars*2));
      assert((unsigned)Yindex<(n_vert_vars*2));
      AddValA(Xindex,Yindex,val[i][j]);
      AddValA(Xindex+1,Yindex+1,val[i][j]);
    }
}

///add a 3x3 block matrix to the system matrix...
///indexes are specified in the 3x3 matrix of x,y pairs
///indexes must be multiplied by 2 cause u and v
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::Add44Block(double val[4][4],int index[4][4][2])
{
  for (int i=0;i<4;i++)
    for (int j=0;j<4;j++)
    {
      ///add for both u and v
      int Xindex=index[i][j][0]*2;
      int Yindex=index[i][j][1]*2;
      assert((unsigned)Xindex<(n_vert_vars*2));
      assert((unsigned)Yindex<(n_vert_vars*2));
      AddValA(Xindex,Yindex,val[i][j]);
      AddValA(Xindex+1,Yindex+1,val[i][j]);
    }
}
///END SYSTEM ACCESS METHODS

///START COMMON MATH FUNCTIONS
///return the complex encoding the rotation
///for a given missmatch interval
template <typename DerivedV, typename DerivedF>
IGL_INLINE std::complex<double> igl::comiso::PoissonSolver<DerivedV, DerivedF>::GetRotationComplex(int interval)
{
  assert((interval>=0)&&(interval<4));

  switch(interval)
  {
    case 0:return std::complex<double>(1,0);
    case 1:return std::complex<double>(0,1);
    case 2:return std::complex<double>(-1,0);
    default:return std::complex<double>(0,-1);
  }
}

///END COMMON MATH FUNCTIONS


///START ENERGY MINIMIZATION PART
///initialize the LHS for a given face
///for minimization of Dirichlet's energy
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::perElementLHS(int f,
                                                                      double val[3][3],
                                                                      int index[3][3][2])
{
  ///initialize to zero
  for (int x=0;x<3;x++)
    for (int y=0;y<3;y++)
      val[x][y]=0;

  ///get the vertices
  int v[3];
  v[0] = F(f,0);
  v[1] = F(f,1);
  v[2] = F(f,2);

  ///get the indexes of vertex instance (to consider cuts)
  ///for the current face
  int Vindexes[3];
  Vindexes[0]=HandleS_Index(f,0);
  Vindexes[1]=HandleS_Index(f,1);
  Vindexes[2]=HandleS_Index(f,2);

  ///initialize the indexes for the block
  for (int x=0;x<3;x++)
    for (int y=0;y<3;y++)
    {
      index[x][y][0]=Vindexes[x];
      index[x][y][1]=Vindexes[y];
    }

  ///initialize edges
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> e[3];
  for (int k=0;k<3;k++)
    e[k] = V.row(v[(k+2)%3]) - V.row(v[(k+1)%3]);

  ///then consider area but also considering scale factor dur to overlaps

  double areaT = doublearea(f)/2.0;

  for (int x=0;x<3;x++)
    for (int y=0;y<3;y++)
      if (x!=y)
      {
        double num =  (e[x].dot(e[y]));
        val[x][y]  =  num/(4.0*areaT);
        val[x][y]  *= Handle_Stiffness[f];//f->stiffening;
      }

  ///set the matrix as diagonal
  SetDiagonal(val);
}

///initialize the RHS for a given face
///for minimization of Dirichlet's energy
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::perElementRHS(int f,
                                                                      double b[6],
                                                                      double vector_field_scale)
{

  /// then set the rhs
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> scaled_Kreal;
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> scaled_Kimag;
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> fNorm = N.row(f);
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> p[3];
  p[0] = V.row(F(f,0));
  p[1] = V.row(F(f,1));
  p[2] = V.row(F(f,2));

  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> neg_t[3];
  neg_t[0] = fNorm.cross(p[2] - p[1]);
  neg_t[1] = fNorm.cross(p[0] - p[2]);
  neg_t[2] = fNorm.cross(p[1] - p[0]);

  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> K1,K2;
  K1 = PD1.row(f);
  K2 = -PD2.row(f); // TODO: the "-" accounts for the orientation of local_basis.h, adapt the code before and remove the "-"

  scaled_Kreal = K1*(vector_field_scale)/2;
  scaled_Kimag = K2*(vector_field_scale)/2;

  double stiff_val = Handle_Stiffness[f];

  b[0] = scaled_Kreal.dot(neg_t[0]) * stiff_val;
  b[1] = scaled_Kimag.dot(neg_t[0]) * stiff_val;
  b[2] = scaled_Kreal.dot(neg_t[1]) * stiff_val;
  b[3] = scaled_Kimag.dot(neg_t[1]) * stiff_val;
  b[4] = scaled_Kreal.dot(neg_t[2]) * stiff_val;
  b[5] = scaled_Kimag.dot(neg_t[2]) * stiff_val;

  //    if (f == 0)
  //    {
  //      cerr << "DEBUG!!!" << endl;
  //
  //
  //      for (unsigned z = 0; z<6; ++z)
  //        cerr << b[z] << " ";
  //      cerr << endl;
  //
  //      scaled_Kreal = K1*(vector_field_scale)/2;
  //      scaled_Kimag = -K2*(vector_field_scale)/2;
  //
  //      double stiff_val = Handle_Stiffness[f];
  //
  //      b[0] = scaled_Kreal.dot(neg_t[0]) * stiff_val;
  //      b[1] = scaled_Kimag.dot(neg_t[0]) * stiff_val;
  //      b[2] = scaled_Kreal.dot(neg_t[1]) * stiff_val;
  //      b[3] = scaled_Kimag.dot(neg_t[1]) * stiff_val;
  //      b[4] = scaled_Kreal.dot(neg_t[2]) * stiff_val;
  //      b[5] = scaled_Kimag.dot(neg_t[2]) * stiff_val;
  //
  //      for (unsigned z = 0; z<6; ++z)
  //        cerr << b[z] << " ";
  //      cerr << endl;
  //
  //    }

}

///evaluate the LHS and RHS for a single face
///for minimization of Dirichlet's energy
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::PerElementSystemReal(int f,
                                                                             double val[3][3],
                                                                             int index[3][3][2],
                                                                             double b[6],
                                                                             double vector_field_scale)
{
  perElementLHS(f,val,index);
  perElementRHS(f,b,vector_field_scale);
}
///END ENERGY MINIMIZATION PART

///START FIXING VERTICES
///set a given vertex as fixed
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddFixedVertex(int v)
{
  n_fixed_vars++;
  Hard_constraints.push_back(v);
}

///find vertex to fix in case we're using
///a vector field NB: multiple components not handled
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::FindFixedVertField()
{
  Hard_constraints.clear();

  n_fixed_vars=0;
  //fix the first singularity
  for (unsigned int v=0;v<V.rows();v++)
  {
    if (Handle_Singular(v))
    {
      AddFixedVertex(v);
      UV.row(v) << 0,0;
      return;
    }
  }

  ///if anything fixed fix the first
  AddFixedVertex(0); // TODO HERE IT ISSSSSS
  UV.row(0) << 0,0;
  std::cerr << "No vertices to fix, I am fixing the first vertex to the origin!" << std::endl;
}

///find hard constraint depending if using or not
///a vector field
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::FindFixedVert()
{
  Hard_constraints.clear();
  FindFixedVertField();
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE int igl::comiso::PoissonSolver<DerivedV, DerivedF>::GetFirstVertexIndex(int v)
{
  return HandleS_Index(VF[v][0],VFi[v][0]);
}

///fix the vertices which are flagged as fixed
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::FixBlockedVertex()
{
  int offset_row = n_vert_vars*2 + num_cut_constraint*2;

  unsigned int constr_num = 0;
  for (unsigned int i=0;i<Hard_constraints.size();i++)
  {
    int v = Hard_constraints[i];

    ///get first index of the vertex that must blocked
    //int index=v->vertex_index[0];
    int index = GetFirstVertexIndex(v);

    ///multiply times 2 because of uv
    int indexvert = index*2;

    ///find the first free row to add the constraint
    int indexRow = (offset_row+constr_num*2);
    int indexCol = indexRow;

    ///add fixing constraint LHS
    AddValA(indexRow,indexvert,1);
    AddValA(indexRow+1,indexvert+1,1);

    ///add fixing constraint RHS
    AddValB(indexCol,  UV(v,0));
    AddValB(indexCol+1,UV(v,1));

    constr_num++;
  }
  assert(constr_num==n_fixed_vars);
}
///END FIXING VERTICES

///HANDLING SINGULARITY
//set the singularity round to integer location
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddSingularityRound()
{
  for (unsigned int v=0;v<V.rows();v++)
  {
    if (Handle_Singular(v))
    {
      int index0=GetFirstVertexIndex(v);
      ids_to_round.push_back( index0*2   );
      ids_to_round.push_back((index0*2)+1);
    }
  }
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AddToRoundVertices(std::vector<int> ids)
{
  for (size_t i = 0; i < ids.size(); ++i)
  {
    if (ids[i] < 0 || ids[i] >= V.rows())
      std::cerr << "WARNING: Ignored round vertex constraint, vertex " << ids[i] << " does not exist in the mesh." << std::endl;
    int index0 = GetFirstVertexIndex(ids[i]);
    ids_to_round.push_back( index0*2   );
    ids_to_round.push_back((index0*2)+1);
  }
}

///START GENERIC SYSTEM FUNCTIONS
//build the laplacian matrix cyclyng over all rangemaps
//and over all faces
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::BuildLaplacianMatrix(double vfscale)
{
  ///then for each face
  for (unsigned int f=0;f<F.rows();f++)
  {

    int var_idx[3]; //vertex variable indices

    for(int k = 0; k < 3; ++k)
      var_idx[k] = HandleS_Index(f,k);

    ///block of variables
    double val[3][3];
    ///block of vertex indexes
    int index[3][3][2];
    ///righe hand side
    double b[6];
    ///compute the system for the given face
    PerElementSystemReal(f, val,index, b, vfscale);

    //Add the element to the matrix
    Add33Block(val,index);

    ///add right hand side
    AddRHS(b,var_idx);
  }
}

///find different sized of the system
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::FindSizes()
{
  ///find the vertex that need to be fixed
  FindFixedVert();

  ///REAL PART
  n_vert_vars = Handle_SystemInfo.num_vert_variables;

  ///INTEGER PART
  ///the total number of integer variables
  n_integer_vars = Handle_SystemInfo.num_integer_cuts;

  ///CONSTRAINT PART
  num_cut_constraint = Handle_SystemInfo.EdgeSeamInfo.size()*2;

  num_constraint_equations = num_cut_constraint*2 + n_fixed_vars*2 + num_userdefined_constraint;

  ///total variable of the system
  num_total_vars = n_vert_vars*2+n_integer_vars*2;

  ///initialize matrix size

  system_size = num_total_vars + num_constraint_equations;

  if (DEBUGPRINT)     printf("\n*** SYSTEM VARIABLES *** \n");
  if (DEBUGPRINT)     printf("* NUM REAL VERTEX VARIABLES %d \n",n_vert_vars);

  if (DEBUGPRINT)     printf("\n*** SINGULARITY *** \n ");
  if (DEBUGPRINT)     printf("* NUM SINGULARITY %d\n",(int)ids_to_round.size()/2);

  if (DEBUGPRINT)     printf("\n*** INTEGER VARIABLES *** \n");
  if (DEBUGPRINT)     printf("* NUM INTEGER VARIABLES %d \n",(int)n_integer_vars);

  if (DEBUGPRINT)     printf("\n*** CONSTRAINTS *** \n ");
  if (DEBUGPRINT)     printf("* NUM FIXED CONSTRAINTS %d\n",n_fixed_vars);
  if (DEBUGPRINT)     printf("* NUM CUTS CONSTRAINTS %d\n",num_cut_constraint);
  if (DEBUGPRINT)     printf("* NUM USER DEFINED CONSTRAINTS %d\n",num_userdefined_constraint);

  if (DEBUGPRINT)     printf("\n*** TOTAL SIZE *** \n");
  if (DEBUGPRINT)     printf("* TOTAL VARIABLE SIZE (WITH INTEGER TRASL) %d \n",num_total_vars);
  if (DEBUGPRINT)     printf("* TOTAL CONSTRAINTS %d \n",num_constraint_equations);
  if (DEBUGPRINT)     printf("* MATRIX SIZE  %d \n",system_size);
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::AllocateSystem()
{
  S.initialize(system_size, system_size);
  printf("\n INITIALIZED SPARSE MATRIX OF %d x %d \n",system_size, system_size);
}

///intitialize the whole matrix
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::InitMatrix()
{
  FindSizes();
  AllocateSystem();
}

///map back coordinates after that
///the system has been solved
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::MapCoords()
{
  ///map coords to faces
  for (unsigned int f=0;f<F.rows();f++)
  {

    for (int k=0;k<3;k++)
    {
      //get the index of the variable in the system
      int indexUV = HandleS_Index(f,k);
      ///then get U and V coords
      double U=X[indexUV*2];
      double V=X[indexUV*2+1];

      WUV(f,k*2 + 0) = U;
      WUV(f,k*2 + 1) = V;
    }
  }

#if 0
  ///initialize the vector of integer variables to return their values
  Handle_SystemInfo.IntegerValues.resize(n_integer_vars*2);
  int baseIndex = (n_vert_vars)*2;
  int endIndex  = baseIndex+n_integer_vars*2;
  int index     = 0;
  for (int i=baseIndex; i<endIndex; i++)
  {
    ///assert that the value is an integer value
    double value=X[i];
    double diff = value-(int)floor(value+0.5);
    assert(diff<0.00000001);
    Handle_SystemInfo.IntegerValues[index] = value;
    index++;
  }
#endif
}

///END GENERIC SYSTEM FUNCTIONS

///set the constraints for the inter-range cuts
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::BuildSeamConstraintsExplicitTranslation()
{
  ///add constraint(s) for every seam edge (not halfedge)
  int offset_row = n_vert_vars;
  ///current constraint row
  int constr_row = offset_row;
  ///current constraint
  unsigned int constr_num = 0;

  for (unsigned int i=0; i<num_cut_constraint/2; i++)
  {
    unsigned char interval = Handle_SystemInfo.EdgeSeamInfo[i].MMatch;
    if (interval==1)
      interval=3;
    else
      if(interval==3)
        interval=1;

    int p0  = Handle_SystemInfo.EdgeSeamInfo[i].v0;
    int p1  = Handle_SystemInfo.EdgeSeamInfo[i].v1;
    int p0p = Handle_SystemInfo.EdgeSeamInfo[i].v0p;
    int p1p = Handle_SystemInfo.EdgeSeamInfo[i].v1p;

    std::complex<double> rot = GetRotationComplex(interval);

    ///get the integer variable
    int integerVar = offset_row+Handle_SystemInfo.EdgeSeamInfo[i].integerVar;

    if (integer_rounding)
    {
      ids_to_round.push_back(integerVar*2);
      ids_to_round.push_back(integerVar*2+1);
    }

    AddComplexA(constr_row, p0 ,  rot);
    AddComplexA(constr_row, p0p,   -1);
    ///then translation...considering the rotation
    ///due to substitution
    AddComplexA(constr_row, integerVar, 1);

    AddValB(2*constr_row  ,0);
    AddValB(2*constr_row+1,0);
    constr_row +=1;
    constr_num++;

    AddComplexA(constr_row, p1,  rot);
    AddComplexA(constr_row, p1p, -1);

    ///other translation
    AddComplexA(constr_row, integerVar  , 1);

    AddValB(2*constr_row,0);
    AddValB(2*constr_row+1,0);

    constr_row +=1;
    constr_num++;
  }
}

///set the constraints for the inter-range cuts
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::BuildUserDefinedConstraints()
{
  /// the user defined constraints are at the end
  int offset_row = n_vert_vars*2 + num_cut_constraint*2 + n_fixed_vars*2;

  ///current constraint row
  int constr_row = offset_row;

  assert(num_userdefined_constraint == userdefined_constraints.size());

  for (unsigned int i=0; i<num_userdefined_constraint; i++)
  {
    for (unsigned int j=0; j<userdefined_constraints[i].size()-1; ++j)
    {
      AddValA(constr_row, j ,  userdefined_constraints[i][j]);
    }

    AddValB(constr_row,userdefined_constraints[i][userdefined_constraints[i].size()-1]);

    constr_row +=1;
  }
}

///call of the mixed integer solver
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::MixedIntegerSolve(double cone_grid_res,
                                                                          bool direct_round,
                                                                          int localIter)
{
  X = std::vector<double>((n_vert_vars+n_integer_vars)*2);

  ///variables part
  int ScalarSize = n_vert_vars*2;
  int SizeMatrix = (n_vert_vars+n_integer_vars)*2;

  if (DEBUGPRINT)
    printf("\n ALLOCATED X \n");

  ///matrix A
  gmm::col_matrix< gmm::wsvector< double > > A(SizeMatrix,SizeMatrix); // lhs matrix variables +

  ///constraints part
  int CsizeX = num_constraint_equations;
  int CsizeY = SizeMatrix+1;
  gmm::row_matrix< gmm::wsvector< double > > C(CsizeX,CsizeY); // constraints

  if (DEBUGPRINT)
    printf("\n ALLOCATED QMM STRUCTURES \n");

  std::vector<double> rhs(SizeMatrix,0);  // rhs

  if (DEBUGPRINT)
    printf("\n ALLOCATED RHS STRUCTURES \n");

  //// copy LHS
  for(int i = 0; i < (int)S.A().nentries(); ++i)
  {
    int row  = S.A().rowind()[i];
    int col  = S.A().colind()[i];
    int size =(int)S.nrows();
    assert(0 <= row && row < size);
    assert(0 <= col && col < size);

    // it's either part of the matrix
    if (row < ScalarSize)
    {
      A(row, col) += S.A().vals()[i];
    }
    // or it's a part of the constraint
    else
    {
      assert ((unsigned int)row < (n_vert_vars+num_constraint_equations)*2);
      int r = row - ScalarSize;
      assert(r   < CsizeX);
      assert(col < CsizeY);
      C(r  , col  ) +=  S.A().vals()[i];
    }
  }

  if (DEBUGPRINT)
    printf("\n SET %d INTEGER VALUES \n",n_integer_vars);

  ///add penalization term for integer variables
  double penalization = 0.000001;
  int offline_index   = ScalarSize;
  for(unsigned int i = 0; i < (n_integer_vars)*2; ++i)
  {
    int index=offline_index+i;
    A(index,index)=penalization;
  }

  if (DEBUGPRINT)
    printf("\n SET RHS \n");

  // copy RHS
  for(int i = 0; i < (int)ScalarSize; ++i)
  {
    rhs[i] = S.getRHSReal(i) * cone_grid_res;
  }

  // copy constraint RHS
  if (DEBUGPRINT)
    printf("\n SET %d CONSTRAINTS \n",num_constraint_equations);

  for(unsigned int i = 0; i < num_constraint_equations; ++i)
  {
    C(i, SizeMatrix) = -S.getRHSReal(ScalarSize + i) * cone_grid_res;
  }

  ///copy values back into S
  COMISO::ConstrainedSolver solver;

  solver.misolver().set_local_iters(localIter);

  solver.misolver().set_direct_rounding(direct_round);

  std::sort(ids_to_round.begin(),ids_to_round.end());
  std::vector<int>::iterator new_end=std::unique(ids_to_round.begin(),ids_to_round.end());
  int dist=distance(ids_to_round.begin(),new_end);
  ids_to_round.resize(dist);

  solver.solve( C, A, X, rhs, ids_to_round, 0.0, false, false);
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::clearUserConstraint()
{
  num_userdefined_constraint = 0;
  userdefined_constraints.clear();
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comiso::PoissonSolver<DerivedV, DerivedF>::addSharpEdgeConstraint(int fid, int vid)
{
  // prepare constraint
  std::vector<int> c(Handle_SystemInfo.num_vert_variables*2 + 1);

  for (size_t i = 0; i < c.size(); ++i)
  {
    c[i] = 0;
  }

  int v1 = F(fid,vid);
  int v2 = F(fid,(vid+1)%3);

  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> e = V.row(v2) - V.row(v1);
  e = e.normalized();

  int v1i = HandleS_Index(fid,vid);//GetFirstVertexIndex(v1);
  int v2i = HandleS_Index(fid,(vid+1)%3);//GetFirstVertexIndex(v2);

  double d1 = fabs(e.dot(PD1.row(fid).normalized()));
  double d2 = fabs(e.dot(PD2.row(fid).normalized()));

  int offset = 0;

  if (d1>d2)
    offset = 1;

  ids_to_round.push_back((v1i * 2) + offset);
  ids_to_round.push_back((v2i * 2) + offset);

  // add constraint
  c[(v1i * 2) + offset] =  1;
  c[(v2i * 2) + offset] = -1;

  // add to the user-defined constraints
  num_userdefined_constraint++;
  userdefined_constraints.push_back(c);

}



template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::MIQ_class(const Eigen::PlainObjectBase<DerivedV> &V_,
                                                                   const Eigen::PlainObjectBase<DerivedF> &F_,
                                                                   const Eigen::PlainObjectBase<DerivedV> &PD1_combed,
                                                                   const Eigen::PlainObjectBase<DerivedV> &PD2_combed,
                                                                   // const Eigen::PlainObjectBase<DerivedV> &BIS1_combed,
                                                                   // const Eigen::PlainObjectBase<DerivedV> &BIS2_combed,
                                                                   const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_MMatch,
                                                                   const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular,
                                                                   // const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_SingularDegree,
                                                                   const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_Seams,
                                                                   Eigen::PlainObjectBase<DerivedU> &UV,
                                                                   Eigen::PlainObjectBase<DerivedF> &FUV,
                                                                   double GradientSize,
                                                                   double Stiffness,
                                                                   bool DirectRound,
                                                                   int iter,
                                                                   int localIter,
                                                                   bool DoRound,
                                                                   bool SingularityRound,
                                                                   std::vector<int> roundVertices,
                                                                   std::vector<std::vector<int> > hardFeatures):
V(V_),
F(F_)
{
  igl::local_basis(V,F,B1,B2,B3);
  igl::triangle_triangle_adjacency(V,F,TT,TTi);

  // Prepare indexing for the linear system
  VertexIndexing<DerivedV, DerivedF> VInd(V, F, TT, TTi, /*BIS1_combed, BIS2_combed,*/ Handle_MMatch, /*Handle_Singular, Handle_SingularDegree,*/ Handle_Seams);

  VInd.InitMapping();
  VInd.InitFaceIntegerVal();
  VInd.InitSeamInfo();

  // Eigen::PlainObjectBase<DerivedV> PD1_combed_for_poisson, PD2_combed_for_poisson;
  // // Rotate by 90 degrees CCW
  // PD1_combed_for_poisson.setZero(PD1_combed.rows(),3);
  // PD2_combed_for_poisson.setZero(PD2_combed.rows(),3);
  // for (unsigned i=0; i<PD1_combed.rows();++i)
  // {
  //   double n1 = PD1_combed.row(i).norm();
  //   double n2 = PD2_combed.row(i).norm();
  //
  //   double a1 = atan2(B2.row(i).dot(PD1_combed.row(i)),B1.row(i).dot(PD1_combed.row(i)));
  //   double a2 = atan2(B2.row(i).dot(PD2_combed.row(i)),B1.row(i).dot(PD2_combed.row(i)));
  //
  //   // a1 += M_PI/2;
  //   // a2 += M_PI/2;
  //
  //
  //   PD1_combed_for_poisson.row(i) = cos(a1) * B1.row(i) + sin(a1) * B2.row(i);
  //   PD2_combed_for_poisson.row(i) = cos(a2) * B1.row(i) + sin(a2) * B2.row(i);
  //
  //   PD1_combed_for_poisson.row(i) = PD1_combed_for_poisson.row(i).normalized() * n1;
  //   PD2_combed_for_poisson.row(i) = PD2_combed_for_poisson.row(i).normalized() * n2;
  // }


  // Assemble the system and solve
  PoissonSolver<DerivedV, DerivedF> PSolver(V,
                                            F,
                                            TT,
                                            TTi,
                                            PD1_combed,
                                            PD2_combed,
                                            VInd.HandleS_Index,
                                            /*VInd.Handle_Singular*/Handle_Singular,
                                            VInd.Handle_SystemInfo);
  Handle_Stiffness = Eigen::VectorXd::Constant(F.rows(),1);


  if (iter > 0) // do stiffening
  {
    for (int i=0;i<iter;i++)
    {
      PSolver.SolvePoisson(Handle_Stiffness, GradientSize,1.f,DirectRound,localIter,DoRound,SingularityRound,roundVertices,hardFeatures);
      int nflips=NumFlips(PSolver.WUV);
      bool folded = updateStiffeningJacobianDistorsion(GradientSize,PSolver.WUV);
      printf("ITERATION %d FLIPS %d \n",i,nflips);
      if (!folded)break;
    }
  }
  else
  {
    PSolver.SolvePoisson(Handle_Stiffness,GradientSize,1.f,DirectRound,localIter,DoRound,SingularityRound,roundVertices,hardFeatures);
  }

  int nflips=NumFlips(PSolver.WUV);
  printf("**** END OPTIMIZING #FLIPS %d  ****\n",nflips);

  fflush(stdout);
  WUV = PSolver.WUV;

}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE void igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::extractUV(Eigen::PlainObjectBase<DerivedU> &UV_out,
                                                                        Eigen::PlainObjectBase<DerivedF> &FUV_out)
{
  //      int f = F.rows();
  int f = WUV.rows();

  unsigned vtfaceid[f*3];
  std::vector<double> vtu;
  std::vector<double> vtv;

  std::vector<std::vector<double> > listUV;
  unsigned counter = 0;

  for (unsigned i=0; i<f; ++i)
  {
    for (unsigned j=0; j<3; ++j)
    {
      std::vector<double> t(3);
      t[0] = WUV(i,j*2 + 0);
      t[1] = WUV(i,j*2 + 1);
      t[2] = counter++;
      listUV.push_back(t);
    }
  }
  std::sort(listUV.begin(),listUV.end());

  counter = 0;
  unsigned k = 0;
  while (k < f*3)
  {
    double u = listUV[k][0];
    double v = listUV[k][1];
    unsigned id = round(listUV[k][2]);

    vtfaceid[id] = counter;
    vtu.push_back(u);
    vtv.push_back(v);

    unsigned j=1;
    while(k+j < f*3 && u == listUV[k+j][0] && v == listUV[k+j][1])
    {
      unsigned tid = round(listUV[k+j][2]);
      vtfaceid[tid] = counter;
      ++j;
    }
    k = k+j;
    counter++;
  }

  UV_out.resize(vtu.size(),2);
  for (unsigned i=0; i<vtu.size(); ++i)
  {
    UV_out(i,0) = vtu[i];
    UV_out(i,1) = vtv[i];
  }

  FUV_out.resize(f,3);

  unsigned vcounter = 0;
  for (unsigned i=0; i<f; ++i)
  {
    FUV_out(i,0)  = vtfaceid[vcounter++];
    FUV_out(i,1)  = vtfaceid[vcounter++];
    FUV_out(i,2)  = vtfaceid[vcounter++];
  }

}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE int igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::NumFlips(const Eigen::MatrixXd& WUV)
{
  int numFl=0;
  for (unsigned int i=0;i<F.rows();i++)
  {
    if (IsFlipped(i, WUV))
      numFl++;
  }
  return numFl;
}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE double igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::Distortion(int f, double h, const Eigen::MatrixXd& WUV)
{
  assert(h > 0);

  Eigen::Vector2d uv0,uv1,uv2;

  uv0 << WUV(f,0), WUV(f,1);
  uv1 << WUV(f,2), WUV(f,3);
  uv2 << WUV(f,4), WUV(f,5);

  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> p0 = V.row(F(f,0));
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> p1 = V.row(F(f,1));
  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> p2 = V.row(F(f,2));

  Eigen::Matrix<typename DerivedV::Scalar, 3, 1> norm = (p1 - p0).cross(p2 - p0);
  double area2 = norm.norm();
  double area2_inv  = 1.0 / area2;
  norm *= area2_inv;

  if (area2 > 0)
  {
    // Singular values of the Jacobian
    Eigen::Matrix<typename DerivedV::Scalar, 3, 1> neg_t0 = norm.cross(p2 - p1);
    Eigen::Matrix<typename DerivedV::Scalar, 3, 1> neg_t1 = norm.cross(p0 - p2);
    Eigen::Matrix<typename DerivedV::Scalar, 3, 1> neg_t2 = norm.cross(p1 - p0);

    Eigen::Matrix<typename DerivedV::Scalar, 3, 1> diffu =  (neg_t0 * uv0(0) +neg_t1 *uv1(0) +  neg_t2 * uv2(0) )*area2_inv;
    Eigen::Matrix<typename DerivedV::Scalar, 3, 1> diffv = (neg_t0 * uv0(1) + neg_t1*uv1(1) +  neg_t2*uv2(1) )*area2_inv;

    // first fundamental form
    double I00 = diffu.dot(diffu);  // guaranteed non-neg
    double I01 = diffu.dot(diffv);  // I01 = I10
    double I11 = diffv.dot(diffv);  // guaranteed non-neg

    // eigenvalues of a 2x2 matrix
    // [a00 a01]
    // [a10 a11]
    // 1/2 * [ (a00 + a11) +/- sqrt((a00 - a11)^2 + 4 a01 a10) ]
    double trI = I00 + I11;                     // guaranteed non-neg
    double diffDiag = I00 - I11;                // guaranteed non-neg
    double sqrtDet = sqrt(std::max(0.0, diffDiag*diffDiag +
                                   4 * I01 * I01)); // guaranteed non-neg
    double sig1 = 0.5 * (trI + sqrtDet); // higher singular value
    double sig2 = 0.5 * (trI - sqrtDet); // lower singular value

    // Avoid sig2 < 0 due to numerical error
    if (fabs(sig2) < 1.0e-8)
      sig2 = 0;

    assert(sig1 >= 0);
    assert(sig2 >= 0);

    if (sig2 < 0) {
      printf("Distortion will be NaN! sig1^2 is negative (%lg)\n",
             sig2);
    }

    // The singular values of the Jacobian are the sqrts of the
    // eigenvalues of the first fundamental form.
    sig1 = sqrt(sig1);
    sig2 = sqrt(sig2);

    // distortion
    double tao = IsFlipped(f,WUV) ? -1 : 1;
    double factor = tao / h;
    double lam = fabs(factor * sig1 - 1) + fabs(factor * sig2 - 1);
    return lam;
  }
  else {
    return 10; // something "large"
  }
}

////////////////////////////////////////////////////////////////////////////
// Approximate the distortion laplacian using a uniform laplacian on
//  the dual mesh:
//      ___________
//      \-1 / \-1 /
//       \ / 3 \ /
//        \-----/
//         \-1 /
//          \ /
//
//  @param[in]  f   facet on which to compute distortion laplacian
//  @param[in]  h   scaling factor applied to cross field
//  @return     distortion laplacian for f
///////////////////////////////////////////////////////////////////////////
template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE double igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::LaplaceDistortion(const int f, double h, const Eigen::MatrixXd& WUV)
{
  double mydist = Distortion(f, h, WUV);
  double lapl=0;
  for (int i=0;i<3;i++)
  {
    if (TT(f,i) != -1)
      lapl += (mydist - Distortion(TT(f,i), h, WUV));
  }
  return lapl;
}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE bool igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::updateStiffeningJacobianDistorsion(double grad_size, const Eigen::MatrixXd& WUV)
{
  bool flipped = NumFlips(WUV)>0;

  if (!flipped)
    return false;

  double maxL=0;
  double maxD=0;

  if (flipped)
  {
    const double c = 1.0;
    const double d = 5.0;

    for (unsigned int i = 0; i < F.rows(); ++i)
    {
      double dist=Distortion(i,grad_size,WUV);
      if (dist > maxD)
        maxD=dist;

      double absLap=fabs(LaplaceDistortion(i, grad_size,WUV));
      if (absLap > maxL)
        maxL = absLap;

      double stiffDelta = std::min(c * absLap, d);

      Handle_Stiffness[i]+=stiffDelta;
    }
  }
  printf("Maximum Distorsion %4.4f \n",maxD);
  printf("Maximum Laplacian %4.4f \n",maxL);
  return flipped;
}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE bool igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::IsFlipped(const Eigen::Vector2d &uv0,
                                                                        const Eigen::Vector2d &uv1,
                                                                        const Eigen::Vector2d &uv2)
{
  Eigen::Vector2d e0 = (uv1-uv0);
  Eigen::Vector2d e1 = (uv2-uv0);

  double Area = e0(0)*e1(1) - e0(1)*e1(0);
  return (Area<=0);
}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE bool igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU>::IsFlipped(
  const int i, const Eigen::MatrixXd& WUV)
{
  Eigen::Vector2d uv0,uv1,uv2;
  uv0 << WUV(i,0), WUV(i,1);
  uv1 << WUV(i,2), WUV(i,3);
  uv2 << WUV(i,4), WUV(i,5);

  return (IsFlipped(uv0,uv1,uv2));
}




template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE void igl::comiso::miq(
  const Eigen::PlainObjectBase<DerivedV> &V,
  const Eigen::PlainObjectBase<DerivedF> &F,
  const Eigen::PlainObjectBase<DerivedV> &PD1_combed,
  const Eigen::PlainObjectBase<DerivedV> &PD2_combed,
  //  const Eigen::PlainObjectBase<DerivedV> &BIS1_combed,
  //  const Eigen::PlainObjectBase<DerivedV> &BIS2_combed,
  const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_MMatch,
  const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular,
  //  const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_SingularDegree,
  const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_Seams,
  Eigen::PlainObjectBase<DerivedU> &UV,
  Eigen::PlainObjectBase<DerivedF> &FUV,
  double GradientSize,
  double Stiffness,
  bool DirectRound,
  int iter,
  int localIter,
  bool DoRound,
  bool SingularityRound,
  std::vector<int> roundVertices,
  std::vector<std::vector<int> > hardFeatures)
{
  GradientSize = GradientSize/(V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();

  igl::comiso::MIQ_class<DerivedV, DerivedF, DerivedU> miq(V,
    F,
    PD1_combed,
    PD2_combed,
    //  BIS1_combed,
    //  BIS2_combed,
    Handle_MMatch,
    Handle_Singular,
    //  Handle_SingularDegree,
    Handle_Seams,
    UV,
    FUV,
    GradientSize,
    Stiffness,
    DirectRound,
    iter,
    localIter,
    DoRound,
    SingularityRound,
    roundVertices,
    hardFeatures);

  miq.extractUV(UV,FUV);
}

template <typename DerivedV, typename DerivedF, typename DerivedU>
IGL_INLINE void igl::comiso::miq(
    const Eigen::PlainObjectBase<DerivedV> &V,
    const Eigen::PlainObjectBase<DerivedF> &F,
    const Eigen::PlainObjectBase<DerivedV> &PD1,
    const Eigen::PlainObjectBase<DerivedV> &PD2,
    Eigen::PlainObjectBase<DerivedU> &UV,
    Eigen::PlainObjectBase<DerivedF> &FUV,
    double GradientSize,
    double Stiffness,
    bool DirectRound,
    int iter,
    int localIter,
    bool DoRound,
    bool SingularityRound,
    std::vector<int> roundVertices,
    std::vector<std::vector<int> > hardFeatures)
{
  // Eigen::MatrixXd PD2i = PD2;
  // if (PD2i.size() == 0)
  // {
  // Eigen::MatrixXd B1, B2, B3;
  // igl::local_basis(V,F,B1,B2,B3);
  // PD2i = igl::rotate_vectors(V,Eigen::VectorXd::Constant(1,M_PI/2),B1,B2);
  // }

  Eigen::PlainObjectBase<DerivedV> BIS1, BIS2;
  igl::compute_frame_field_bisectors(V, F, PD1, PD2, BIS1, BIS2);

  Eigen::PlainObjectBase<DerivedV> BIS1_combed, BIS2_combed;
  igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);

  Eigen::PlainObjectBase<DerivedF> Handle_MMatch;
  igl::cross_field_missmatch(V, F, BIS1_combed, BIS2_combed, true, Handle_MMatch);

  Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;
  igl::find_cross_field_singularities(V, F, Handle_MMatch, isSingularity, singularityIndex);

  Eigen::Matrix<int, Eigen::Dynamic, 3> Handle_Seams;
  igl::cut_mesh_from_singularities(V, F, Handle_MMatch, Handle_Seams);

  Eigen::PlainObjectBase<DerivedV> PD1_combed, PD2_combed;
  igl::comb_frame_field(V, F, PD1, PD2, BIS1_combed, BIS2_combed, PD1_combed, PD2_combed);

  igl::comiso::miq(V,
           F,
           PD1_combed,
           PD2_combed,
           //  BIS1_combed,
           //  BIS2_combed,
           Handle_MMatch,
           isSingularity,
           //  singularityIndex,
           Handle_Seams,
           UV,
           FUV,
           GradientSize,
           Stiffness,
           DirectRound,
           iter,
           localIter,
           DoRound,
           SingularityRound,
           roundVertices,
           hardFeatures);

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::comiso::miq<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, double, double, bool, int, int, bool, bool, std::vector<int, std::allocator<int> >, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >);
template void igl::comiso::miq<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 3, 0, -1, 3> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 3, 0, -1, 3> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, double, double, bool, int, int, bool, bool, std::vector<int, std::allocator<int> >, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >);
template void igl::comiso::miq<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, double, double, bool, int, int, bool, bool, std::vector<int, std::allocator<int> >, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >);
#endif
