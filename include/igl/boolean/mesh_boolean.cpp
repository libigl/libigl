#include "mesh_boolean.h"
#include <igl/peal_outer_hull_layers.h>
#include <igl/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>
#include <igl/unique_simplices.h>
#include <iostream>

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&, 
          Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&)> 
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J);
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC>
IGL_INLINE void igl::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC)
{
  Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1> J;
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&, 
          Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1>&)> 
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J);
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3>&,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&)> 
    & resolve_fun,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MeshBooleanType eff_type = type;
  // Concatenate A and B into a single mesh
  typedef typename DerivedVC::Scalar Scalar;
  typedef typename DerivedFC::Scalar Index;
  typedef Matrix<Scalar,Dynamic,3> MatrixX3S;
  typedef Matrix<Index,Dynamic,3> MatrixX3I;
  typedef Matrix<Index,Dynamic,2> MatrixX2I;
  typedef Matrix<Index,Dynamic,1> VectorXI;
  typedef Matrix<typename DerivedJ::Scalar,Dynamic,1> VectorXJ;
  MatrixX3S V(VA.rows()+VB.rows(),3);
  MatrixX3I F(FA.rows()+FB.rows(),3);
  V.block(0,0,VA.rows(),VA.cols()) = VA;
  V.block(VA.rows(),0,VB.rows(),VB.cols()) = VB;
  switch(type)
  {
    // Minus is implemented by flipping B and computing union
    case MESH_BOOLEAN_TYPE_MINUS:
      F.block(0,0,FA.rows(),FA.cols()) = FA.rowwise().reverse();
      F.block(FA.rows(),0,FB.rows(),FB.cols()) = FB.array()+VA.rows();
      //F.block(0,0,FA.rows(),3) = FA;
      //F.block(FA.rows(),0,FB.rows(),3) = 
      //  FB.rowwise().reverse().array()+VA.rows();
      eff_type = MESH_BOOLEAN_TYPE_INTERSECT;
      break;
    default:
      F.block(0,0,FA.rows(),FA.cols()) = FA;
      F.block(FA.rows(),0,FB.rows(),FB.cols()) = FB.array()+VA.rows();
      break;
  }

  // Resolve intersections (assumes A and B are solid)
  const auto & libigl_resolve = [](
    const MatrixX3S & V,
    const MatrixX3I & F,
    MatrixX3S & CV,
    MatrixX3I & CF,
    VectorXJ & J)
  {
    MatrixX3S SV;
    MatrixX3I SF;
    MatrixX2I SIF;
    VectorXI SIM,UIM;
    RemeshSelfIntersectionsParam params;
    remesh_self_intersections(V,F,params,SV,SF,SIF,J,SIM);
    for_each(SF.data(),SF.data()+SF.size(),[&SIM](int & a){a=SIM(a);});
    {
      remove_unreferenced(SV,SF,CV,CF,UIM);
    }
  };

  MatrixX3S CV;
  MatrixX3I CF;
  VectorXJ CJ;
  if(resolve_fun)
  {
    resolve_fun(V,F,CV,CF,CJ);
  }else
  {
    libigl_resolve(V,F,CV,CF,CJ);
  }

  if(type == MESH_BOOLEAN_TYPE_RESOLVE)
  {
    FC = CF;
    VC = CV;
    J = CJ;
    return;
  }

  Matrix<bool,Dynamic,1> from_A(CF.rows());
  // Peal layers keeping track of odd and even flips
  Matrix<bool,Dynamic,1> odd;
  Matrix<bool,Dynamic,1> flip;
  peal_outer_hull_layers(CV,CF,odd,flip);

  const Index m = CF.rows();
  // Faces of output vG[i] = j means ith face of output should be jth face in F
  std::vector<Index> vG;
  // Whether faces of output should be flipped, Gflip[i] = true means ith face
  // of output should be F.row(vG[i]).reverse() rather than F.row(vG[i])
  std::vector<bool> Gflip;
  for(Index f = 0;f<m;f++)
  {
    switch(eff_type)
    {
      case MESH_BOOLEAN_TYPE_XOR:
      case MESH_BOOLEAN_TYPE_UNION:
        if((odd(f)&&!flip(f))||(!odd(f)&&flip(f)))
        {
          vG.push_back(f);
          Gflip.push_back(false);
        }else if(eff_type == MESH_BOOLEAN_TYPE_XOR)
        {
          vG.push_back(f);
          Gflip.push_back(true);
        }
        break;
      case MESH_BOOLEAN_TYPE_INTERSECT:
        if((!odd(f) && !flip(f)) || (odd(f) && flip(f)))
        {
          vG.push_back(f);
          Gflip.push_back(type == MESH_BOOLEAN_TYPE_MINUS);
        }
        break;
      default:
        assert(false && "Unknown type");
        return;
    }
  }
  const Index gm = vG.size();
  MatrixX3I G(gm,3);
  VectorXi GJ(gm,1);
  for(Index g = 0;g<gm;g++)
  {
    G.row(g) = Gflip[g] ? CF.row(vG[g]).reverse().eval() : CF.row(vG[g]);
    GJ(g) = CJ(vG[g]);
  }
  // remove duplicates: cancel out in all cases? assertion in others?
  {
    MatrixXi oldG = G;
    VectorXi IA,_1;
    unique_simplices(oldG,G,IA,_1);
    assert(IA.rows() == G.rows());
    J.resize(IA.rows(),1);
    for(size_t j = 0;j<(size_t)J.size();j++)
    {
      J(j) = GJ(IA(j));
    }
  }
  // remove unreferenced vertices
  VectorXi newIM;
  remove_unreferenced(CV,G,VC,FC,newIM);
  //cerr<<"warning not removing unref"<<endl;
  //VC = CV;
  //FC = G;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
