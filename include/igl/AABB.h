// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_AABB_H
#define IGL_AABB_H

// Implementation of semi-general purpose axis-aligned bounding box hierarchy.
// The mesh (V,Ele) is stored and managed by the caller and each routine here
// simply takes it as references (it better not change between calls).
//
// It's a little annoying that the Dimension is a template parameter and not
// picked up at run time from V. This leads to duplicated code for 2d/3d (up to
// dim).
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
namespace igl
{
  template <typename DerivedV, int DIM>
    class AABB 
    {
public:
      typedef typename DerivedV::Scalar Scalar;
      typedef Eigen::Matrix<Scalar,1,DIM> RowVectorDIMS;
      typedef Eigen::Matrix<Scalar,DIM,1> VectorDIMS;
      typedef Eigen::Matrix<Scalar,Eigen::Dynamic,DIM> MatrixXDIMS;
      // Shared pointers are slower...
      AABB * m_left;
      AABB * m_right;
      Eigen::AlignedBox<Scalar,DIM> m_box;
      // -1 non-leaf
      int m_primitive;
      //Scalar m_max_sqr_d;
      //int m_depth;
      AABB():
        m_left(NULL), m_right(NULL),
        m_box(), m_primitive(-1)
        //m_max_sqr_d(std::numeric_limits<double>::infinity()),
        //m_depth(0)
    {}
      // http://stackoverflow.com/a/3279550/148668
      AABB(const AABB& other):
        m_left(other.m_left ? new AABB(*other.m_left) : NULL),
        m_right(other.m_right ? new AABB(*other.m_right) : NULL),
        m_box(other.m_box),
        m_primitive(other.m_primitive)
        //m_max_sqr_d(other.m_max_sqr_d),
        //m_depth(std::max(
        //   m_left ? m_left->m_depth + 1 : 0,
        //   m_right ? m_right->m_depth + 1 : 0))
        {
        }
      // copy-swap idiom
      friend void swap(AABB& first, AABB& second)
      {
        // Enable ADL
        using std::swap;
        swap(first.m_left,second.m_left);
        swap(first.m_right,second.m_right);
        swap(first.m_box,second.m_box);
        swap(first.m_primitive,second.m_primitive);
        //swap(first.m_max_sqr_d,second.m_max_sqr_d);
        //swap(first.m_depth,second.m_depth);
      }
      // Pass-by-value (aka copy)
      AABB& operator=(AABB other)
      {
        swap(*this,other);
        return *this;
      }
      AABB(AABB&& other):
        // initialize via default constructor
        AABB() 
    {
      swap(*this,other);
    }
      // Seems like there should have been an elegant solution to this using
      // the copy-swap idiom above:
      inline void deinit()
      {
        m_primitive = -1;
        m_box = Eigen::AlignedBox<Scalar,DIM>();
        delete m_left;
        m_left = NULL;
        delete m_right;
        m_right = NULL;
      }
      ~AABB()
      {
        deinit();
      }
      // Build an Axis-Aligned Bounding Box tree for a given mesh and given
      // serialization of a previous AABB tree.
      //
      // Inputs:
      //   V  #V by dim list of mesh vertex positions. 
      //   Ele  #Ele by dim+1 list of mesh indices into #V. 
      //   bb_mins  max_tree by dim list of bounding box min corner positions
      //   bb_maxs  max_tree by dim list of bounding box max corner positions
      //   elements  max_tree list of element or (not leaf id) indices into Ele
      //   i  recursive call index {0}
      template <typename Derivedbb_mins, typename Derivedbb_maxs>
        inline void init(
            const Eigen::PlainObjectBase<DerivedV> & V,
            const Eigen::MatrixXi & Ele, 
            const Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
            const Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
            const Eigen::VectorXi & elements,
            const int i = 0);
      // Wrapper for root with empty serialization
      inline void init(
          const Eigen::PlainObjectBase<DerivedV> & V,
          const Eigen::MatrixXi & Ele);
      // Build an Axis-Aligned Bounding Box tree for a given mesh.
      //
      // Inputs:
      //   V  #V by dim list of mesh vertex positions. 
      //   Ele  #Ele by dim+1 list of mesh indices into #V. 
      //   SI  #Ele by dim list revealing for each coordinate where Ele's
      //     barycenters would be sorted: SI(e,d) = i --> the dth coordinate of
      //     the barycenter of the eth element would be placed at position i in a
      //     sorted list.
      //   I  #I list of indices into Ele of elements to include (for recursive
      //     calls)
      // 
      inline void init(
          const Eigen::PlainObjectBase<DerivedV> & V,
          const Eigen::MatrixXi & Ele, 
          const Eigen::MatrixXi & SI,
          const Eigen::VectorXi & I);
      // Return whether at leaf node
      inline bool is_leaf() const;
      // Find the indices of elements containing given point.
      //
      // Inputs:
      //   V  #V by dim list of mesh vertex positions. **Should be same as used to
      //     construct mesh.**
      //   Ele  #Ele by dim+1 list of mesh indices into #V. **Should be same as used to
      //     construct mesh.**
      //   q  dim row-vector query position
      //   first  whether to only return first element containing q
      // Returns:
      //   list of indices of elements containing q
      template <typename Derivedq>
      inline std::vector<int> find(
          const Eigen::PlainObjectBase<DerivedV> & V,
          const Eigen::MatrixXi & Ele, 
          const Eigen::PlainObjectBase<Derivedq> & q,
          const bool first=false) const;

      // If number of elements m then total tree size should be 2*h where h is
      // the deepest depth 2^ceil(log(#Ele*2-1))
      inline int subtree_size() const;

      // Serialize this class into 3 arrays (so we can pass it pack to matlab)
      //
      // Outputs:
      //   bb_mins  max_tree by dim list of bounding box min corner positions
      //   bb_maxs  max_tree by dim list of bounding box max corner positions
      //   elements  max_tree list of element or (not leaf id) indices into Ele
      //   i  recursive call index into these arrays {0}
      template <typename Derivedbb_mins, typename Derivedbb_maxs>
        inline void serialize(
            Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
            Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
            Eigen::VectorXi & elements,
            const int i = 0) const;
      // Compute squared distance to a query point
      //
      // Inputs:
      //   V  #V by dim list of vertex positions
      //   Ele  #Ele by dim list of simplex indices
      //   P  3 list of query point coordinates
      //   min_sqr_d  current minimum squared distance (only find distances
      //   less than this)
      // Outputs:
      //   I  #P list of facet indices corresponding to smallest distances
      //   C  #P by 3 list of closest points
      // Returns squared distance
      //
      // Known bugs: currently assumes Elements are triangles regardless of
      // dimension.
      inline Scalar squared_distance(
          const Eigen::PlainObjectBase<DerivedV> & V,
          const Eigen::MatrixXi & Ele, 
          const RowVectorDIMS & p,
          int & i,
          RowVectorDIMS & c) const;
private:
      inline Scalar squared_distance(
          const Eigen::PlainObjectBase<DerivedV> & V,
          const Eigen::MatrixXi & Ele, 
          const RowVectorDIMS & p,
          const Scalar min_sqr_d,
          int & i,
          RowVectorDIMS & c) const;
public:
      template <
        typename DerivedP, 
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      inline void squared_distance(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::MatrixXi & Ele, 
        const Eigen::PlainObjectBase<DerivedP> & P,
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;

      template < 
        typename Derivedother_V,
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      inline void squared_distance(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::MatrixXi & Ele, 
        const AABB<Derivedother_V,DIM> & other,
        const Eigen::PlainObjectBase<Derivedother_V> & other_V,
        const Eigen::MatrixXi & other_Ele, 
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;
private:
      template < 
        typename Derivedother_V,
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      inline Scalar squared_distance_helper(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::MatrixXi & Ele, 
        const AABB<Derivedother_V,DIM> * other,
        const Eigen::PlainObjectBase<Derivedother_V> & other_V,
        const Eigen::MatrixXi & other_Ele, 
        const Scalar min_sqr_d,
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;
      // Helper function for leaves: works in-place on sqr_d
      inline void leaf_squared_distance(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::MatrixXi & Ele, 
        const RowVectorDIMS & p,
        Scalar & sqr_d,
        int & i,
        RowVectorDIMS & c) const;
      inline void set_min(
        const RowVectorDIMS & p,
        const Scalar sqr_d_candidate,
        const int i_candidate,
        const RowVectorDIMS & c_candidate,
        Scalar & sqr_d,
        int & i,
        RowVectorDIMS & c) const;
public:
      static
      inline void barycentric_coordinates(
        const RowVectorDIMS & p, 
        const RowVectorDIMS & a, 
        const RowVectorDIMS & b, 
        const RowVectorDIMS & c,
        Eigen::Matrix<Scalar,1,3> & bary);
public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
}

// Implementation
#include "EPS.h"
#include "barycenter.h"
#include "colon.h"
#include "colon.h"
#include "doublearea.h"
#include "matlab_format.h"
#include "project_to_line_segment.h"
#include "sort.h"
#include "volume.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <list>

template <typename DerivedV, int DIM>
  template <typename Derivedbb_mins, typename Derivedbb_maxs>
inline void igl::AABB<DerivedV,DIM>::init(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
    const Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
    const Eigen::VectorXi & elements,
    const int i)
{
  using namespace std;
  using namespace Eigen;
  if(bb_mins.size() > 0)
  {
    assert(bb_mins.rows() == bb_maxs.rows() && "Serial tree arrays must match");
    assert(bb_mins.cols() == V.cols() && "Serial tree array dim must match V");
    assert(bb_mins.cols() == bb_maxs.cols() && "Serial tree arrays must match");
    assert(bb_mins.rows() == elements.rows() &&
        "Serial tree arrays must match");
    // construct from serialization
    m_box.extend(bb_mins.row(i).transpose());
    m_box.extend(bb_maxs.row(i).transpose());
    m_primitive = elements(i);
    // Not leaf then recurse
    if(m_primitive == -1)
    {
      m_left = new AABB();
      m_left->init( V,Ele,bb_mins,bb_maxs,elements,2*i+1);
      m_right = new AABB();
      m_right->init( V,Ele,bb_mins,bb_maxs,elements,2*i+2);
      //m_depth = std::max( m_left->m_depth, m_right->m_depth)+1;
    }
  }else
  {
    VectorXi allI = colon<int>(0,Ele.rows()-1);
    MatrixXDIMS BC;
    if(Ele.cols() == 1)
    {
      // points
      BC = V;
    }else
    {
      // Simplices
      barycenter(V,Ele,BC);
    }
    MatrixXi SI(BC.rows(),BC.cols());
    {
      MatrixXDIMS _;
      MatrixXi IS;
      igl::sort(BC,1,true,_,IS);
      // Need SI(i) to tell which place i would be sorted into
      const int dim = IS.cols();
      for(int i = 0;i<IS.rows();i++)
      {
        for(int d = 0;d<dim;d++)
        {
          SI(IS(i,d),d) = i;
        }
      }
    }
    init(V,Ele,SI,allI);
  }
}

  template <typename DerivedV, int DIM>
inline void igl::AABB<DerivedV,DIM>::init(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::MatrixXi & Ele)
{
  using namespace Eigen;
  return init(V,Ele,MatrixXDIMS(),MatrixXDIMS(),VectorXi(),0);
}

  template <typename DerivedV, int DIM>
inline void igl::AABB<DerivedV,DIM>::init(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::MatrixXi & SI,
    const Eigen::VectorXi & I)
{
  using namespace Eigen;
  using namespace std;
  assert(DIM == V.cols() && "V.cols() should matched declared dimension");
  const Scalar inf = numeric_limits<Scalar>::infinity();
  m_box = AlignedBox<Scalar,DIM>();
  // Compute bounding box
  for(int i = 0;i<I.rows();i++)
  {
    for(int c = 0;c<Ele.cols();c++)
    {
      m_box.extend(V.row(Ele(I(i),c)).transpose());
      m_box.extend(V.row(Ele(I(i),c)).transpose());
    }
  }
  switch(I.size())
  {
    case 0:
      {
        assert(false);
      }
    case 1:
      {
        m_primitive = I(0);
        break;
      }
    default:
      {
        // Compute longest direction
        int max_d = -1;
        m_box.diagonal().maxCoeff(&max_d);
        // Can't use median on BC directly because many may have same value,
        // but can use median on sorted BC indices
        VectorXi SIdI(I.rows());
        for(int i = 0;i<I.rows();i++)
        {
          SIdI(i) = SI(I(i),max_d);
        }
        // Since later I use <= I think I don't need to worry about odd/even
        // Pass by copy to avoid changing input
        const auto median = [](VectorXi A)->Scalar
        {
          size_t n = A.size()/2;
          nth_element(A.data(),A.data()+n,A.data()+A.size());
          if(A.rows() % 2 == 1)
          {
            return A(n);
          }else
          {
            nth_element(A.data(),A.data()+n-1,A.data()+A.size());
            return 0.5*(A(n)+A(n-1));
          }
        };
        const Scalar med = median(SIdI);
        VectorXi LI((I.rows()+1)/2),RI(I.rows()/2);
        assert(LI.rows()+RI.rows() == I.rows());
        // Distribute left and right
        {
          int li = 0;
          int ri = 0;
          for(int i = 0;i<I.rows();i++)
          {
            if(SIdI(i)<=med)
            {
              LI(li++) = I(i);
            }else
            {
              RI(ri++) = I(i);
            }
          }
        }
        //m_depth = 0;
        if(LI.rows()>0)
        {
          m_left = new AABB();
          m_left->init(V,Ele,SI,LI);
          //m_depth = std::max(m_depth, m_left->m_depth+1);
        }
        if(RI.rows()>0)
        {
          m_right = new AABB();
          m_right->init(V,Ele,SI,RI);
          //m_depth = std::max(m_depth, m_right->m_depth+1);
        }
      }
  }
}

template <typename DerivedV, int DIM>
inline bool igl::AABB<DerivedV,DIM>::is_leaf() const
{
  return m_primitive != -1;
}

template <typename DerivedV, int DIM>
template <typename Derivedq>
inline std::vector<int> igl::AABB<DerivedV,DIM>::find(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::PlainObjectBase<Derivedq> & q,
    const bool first) const
{
  using namespace std;
  using namespace Eigen;
  assert(q.size() == DIM && 
      "Query dimension should match aabb dimension");
  assert(Ele.cols() == V.cols()+1 && 
      "AABB::find only makes sense for (d+1)-simplices");
  const Scalar epsilon = igl::EPS<Scalar>();
  // Check if outside bounding box
  bool inside = m_box.contains(q.transpose());
  if(!inside)
  {
    return std::vector<int>();
  }
  assert(m_primitive==-1 || (m_left == NULL && m_right == NULL));
  if(is_leaf())
  {
    // Initialize to some value > -epsilon
    Scalar a1=1,a2=1,a3=1,a4=1;
    switch(DIM)
    {
      case 3:
        {
          // Barycentric coordinates
          typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
          const RowVector3S V1 = V.row(Ele(m_primitive,0));
          const RowVector3S V2 = V.row(Ele(m_primitive,1));
          const RowVector3S V3 = V.row(Ele(m_primitive,2));
          const RowVector3S V4 = V.row(Ele(m_primitive,3));
          a1 = volume_single(V2,V4,V3,(RowVector3S)q);
          a2 = volume_single(V1,V3,V4,(RowVector3S)q);
          a3 = volume_single(V1,V4,V2,(RowVector3S)q);
          a4 = volume_single(V1,V2,V3,(RowVector3S)q);
          break;
        }
      case 2:
        {
          // Barycentric coordinates
          typedef Eigen::Matrix<Scalar,2,1> Vector2S;
          const Vector2S V1 = V.row(Ele(m_primitive,0));
          const Vector2S V2 = V.row(Ele(m_primitive,1));
          const Vector2S V3 = V.row(Ele(m_primitive,2));
          // Hack for now to keep templates simple. If becomes bottleneck
          // consider using std::enable_if_t 
          const Vector2S q2 = q.head(2);
          Scalar a0 = doublearea_single(V1,V2,V3);
          a1 = doublearea_single(V1,V2,q2);
          a2 = doublearea_single(V2,V3,q2);
          a3 = doublearea_single(V3,V1,q2);
          break;
        }
      default:assert(false);
    }
    if(
        a1>=-epsilon && 
        a2>=-epsilon && 
        a3>=-epsilon && 
        a4>=-epsilon)
    {
      return std::vector<int>(1,m_primitive);
    }else
    {
      return std::vector<int>();
    }
  }
  std::vector<int> left = m_left->find(V,Ele,q,first);
  if(first && !left.empty())
  {
    return left;
  }
  std::vector<int> right = m_right->find(V,Ele,q,first);
  if(first)
  {
    return right;
  }
  left.insert(left.end(),right.begin(),right.end());
  return left;
}

template <typename DerivedV, int DIM>
inline int igl::AABB<DerivedV,DIM>::subtree_size() const
{
  // 1 for self
  int n = 1;
  int n_left = 0,n_right = 0;
  if(m_left != NULL)
  {
    n_left = m_left->subtree_size();
  }
  if(m_right != NULL)
  {
    n_right = m_right->subtree_size();
  }
  n += 2*std::max(n_left,n_right);
  return n;
}


template <typename DerivedV, int DIM>
template <typename Derivedbb_mins, typename Derivedbb_maxs>
inline void igl::AABB<DerivedV,DIM>::serialize(
    Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
    Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
    Eigen::VectorXi & elements,
    const int i) const
{
  using namespace std;
  using namespace Eigen;
  // Calling for root then resize output
  if(i==0)
  {
    const int m = subtree_size();
    //cout<<"m: "<<m<<endl;
    bb_mins.resize(m,DIM);
    bb_maxs.resize(m,DIM);
    elements.resize(m,1);
  }
  //cout<<i<<" ";
  bb_mins.row(i) = m_box.min();
  bb_maxs.row(i) = m_box.max();
  elements(i) = m_primitive;
  if(m_left != NULL)
  {
    m_left->serialize(bb_mins,bb_maxs,elements,2*i+1);
  }
  if(m_right != NULL)
  {
    m_right->serialize(bb_mins,bb_maxs,elements,2*i+2);
  }
}

template <typename DerivedV, int DIM>
inline typename igl::AABB<DerivedV,DIM>::Scalar 
igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const RowVectorDIMS & p,
  int & i,
  RowVectorDIMS & c) const
{
  return squared_distance(V,Ele,p,std::numeric_limits<Scalar>::infinity(),i,c);
}


template <typename DerivedV, int DIM>
inline typename igl::AABB<DerivedV,DIM>::Scalar 
igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const RowVectorDIMS & p,
  Scalar min_sqr_d,
  int & i,
  RowVectorDIMS & c) const
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  Scalar sqr_d = min_sqr_d;
  assert(DIM == 3 && "Code has only been tested for DIM == 3");
  assert((Ele.cols() == 3 || Ele.cols() == 2 || Ele.cols() == 1)
    && "Code has only been tested for simplex sizes 3,2,1");

  assert(m_primitive==-1 || (m_left == NULL && m_right == NULL));
  if(is_leaf())
  {
    leaf_squared_distance(V,Ele,p,sqr_d,i,c);
  }else
  {
    bool looked_left = false;
    bool looked_right = false;
    const auto & look_left = [&]()
    {
      int i_left;
      RowVectorDIMS c_left = c;
      Scalar sqr_d_left = m_left->squared_distance(V,Ele,p,sqr_d,i_left,c_left);
      set_min(p,sqr_d_left,i_left,c_left,sqr_d,i,c);
      looked_left = true;
    };
    const auto & look_right = [&]()
    {
      int i_right;
      RowVectorDIMS c_right = c;
      Scalar sqr_d_right = m_right->squared_distance(V,Ele,p,sqr_d,i_right,c_right);
      set_min(p,sqr_d_right,i_right,c_right,sqr_d,i,c);
      looked_right = true;
    };

    // must look left or right if in box
    if(m_left->m_box.contains(p.transpose()))
    {
      look_left();
    }
    if(m_right->m_box.contains(p.transpose()))
    {
      look_right();
    }
    // if haven't looked left and could be less than current min, then look
    Scalar  left_min_sqr_d = m_left->m_box.squaredExteriorDistance(p.transpose());
    Scalar right_min_sqr_d = m_right->m_box.squaredExteriorDistance(p.transpose());
    if(left_min_sqr_d < right_min_sqr_d)
    {
      if(!looked_left && left_min_sqr_d<sqr_d)
      {
        look_left();
      }
      if( !looked_right && right_min_sqr_d<sqr_d)
      {
        look_right();
      }
    }else
    {
      if( !looked_right && right_min_sqr_d<sqr_d)
      {
        look_right();
      }
      if(!looked_left && left_min_sqr_d<sqr_d)
      {
        look_left();
      }
    }
  }
  return sqr_d;
}

template <typename DerivedV, int DIM>
template <
  typename DerivedP, 
  typename DerivedsqrD, 
  typename DerivedI, 
  typename DerivedC>
inline void igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const Eigen::PlainObjectBase<DerivedP> & P,
  Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  assert(P.cols() == V.cols() && "cols in P should match dim of cols in V");
  sqrD.resize(P.rows(),1);
  I.resize(P.rows(),1);
  C.resize(P.rows(),P.cols());
  for(int p = 0;p<P.rows();p++)
  {
    RowVectorDIMS Pp = P.row(p), c;
    int Ip;
    sqrD(p) = squared_distance(V,Ele,Pp,Ip,c);
    I(p) = Ip;
    C.row(p) = c;
  }
}

template <typename DerivedV, int DIM>
template < 
  typename Derivedother_V,
  typename DerivedsqrD, 
  typename DerivedI, 
  typename DerivedC>
inline void igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const AABB<Derivedother_V,DIM> & other,
  const Eigen::PlainObjectBase<Derivedother_V> & other_V,
  const Eigen::MatrixXi & other_Ele, 
  Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  assert(other_Ele.cols() == 1 && 
    "Only implemented for other as list of points");
  assert(other_V.cols() == V.cols() && "other must match this dimension");
  sqrD.setConstant(other_Ele.rows(),1,std::numeric_limits<double>::infinity());
  I.resize(other_Ele.rows(),1);
  C.resize(other_Ele.rows(),other_V.cols());
  // All points in other_V currently think they need to check against root of
  // this. The point of using another AABB is to quickly prune chunks of
  // other_V so that most points just check some subtree of this.

  // This holds a conservative estimate of max(sqr_D) where sqr_D is the
  // current best minimum squared distance for all points in this subtree
  double min_sqr_d = std::numeric_limits<double>::infinity();
  squared_distance_helper(
    V,Ele,&other,other_V,other_Ele,min_sqr_d,sqrD,I,C);
}

template <typename DerivedV, int DIM>
template < 
  typename Derivedother_V,
  typename DerivedsqrD, 
  typename DerivedI, 
  typename DerivedC>
inline typename igl::AABB<DerivedV,DIM>::Scalar igl::AABB<DerivedV,DIM>::squared_distance_helper(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const AABB<Derivedother_V,DIM> * other,
  const Eigen::PlainObjectBase<Derivedother_V> & other_V,
  const Eigen::MatrixXi & other_Ele, 
  const Scalar min_sqr_d,
  Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  using namespace std;
  using namespace Eigen;

  // This implementation is a bit disappointing. There's no major speed up. Any
  // performance gains seem to come from accidental cache coherency and
  // diminish for larger "other" (the opposite of what was intended).

  // Base case
  if(other->is_leaf() && this->is_leaf())
  {
    Scalar sqr_d = sqrD(other->m_primitive);
    int i = I(other->m_primitive);
    RowVectorDIMS c = C.row(      other->m_primitive);
    RowVectorDIMS p = other_V.row(other->m_primitive);
    leaf_squared_distance(V,Ele,p,sqr_d,i,c);
    sqrD( other->m_primitive) = sqr_d;
    I(    other->m_primitive) = i;
    C.row(other->m_primitive) = c;
    //cout<<"leaf: "<<sqr_d<<endl;
    //other->m_max_sqr_d = sqr_d;
    return sqr_d;
  }

  if(other->is_leaf())
  {
    Scalar sqr_d = sqrD(other->m_primitive);
    int i = I(other->m_primitive);
    RowVectorDIMS c = C.row(      other->m_primitive);
    RowVectorDIMS p = other_V.row(other->m_primitive);
    sqr_d = squared_distance(V,Ele,p,sqr_d,i,c);
    sqrD( other->m_primitive) = sqr_d;
    I(    other->m_primitive) = i;
    C.row(other->m_primitive) = c;
    //other->m_max_sqr_d = sqr_d;
    return sqr_d;
  }

  //// Exact minimum squared distance between arbitary primitives inside this and
  //// othre's bounding boxes
  //const auto & min_squared_distance = [&](
  //  const AABB<DerivedV,DIM> * A,
  //  const AABB<Derivedother_V,DIM> * B)->Scalar
  //{
  //  return A->m_box.squaredExteriorDistance(B->m_box);
  //};

  if(this->is_leaf())
  {
    //if(min_squared_distance(this,other) < other->m_max_sqr_d)
    if(true)
    {
      this->squared_distance_helper(
        V,Ele,other->m_left,other_V,other_Ele,0,sqrD,I,C);
      this->squared_distance_helper(
        V,Ele,other->m_right,other_V,other_Ele,0,sqrD,I,C);
    }else
    {
      // This is never reached...
    }
    //// we know other is not a leaf
    //other->m_max_sqr_d = std::max(other->m_left->m_max_sqr_d,other->m_right->m_max_sqr_d);
    return 0;
  }

  // FORCE DOWN TO OTHER LEAF EVAL
  //if(min_squared_distance(this,other) < other->m_max_sqr_d)
  if(true)
  {
    if(true)
    {
      this->squared_distance_helper(
        V,Ele,other->m_left,other_V,other_Ele,0,sqrD,I,C);
      this->squared_distance_helper(
        V,Ele,other->m_right,other_V,other_Ele,0,sqrD,I,C);
    }else // this direction never seems to be faster
    {
      this->m_left->squared_distance_helper(
        V,Ele,other,other_V,other_Ele,0,sqrD,I,C);
      this->m_right->squared_distance_helper(
        V,Ele,other,other_V,other_Ele,0,sqrD,I,C);
    }
  }else
  {
    // this is never reached ... :-(
  }
  //// we know other is not a leaf
  //other->m_max_sqr_d = std::max(other->m_left->m_max_sqr_d,other->m_right->m_max_sqr_d);

  return 0;
#if 0 // False

  // _Very_ conservative approximation of maximum squared distance between
  // primitives inside this and other's bounding boxes
  const auto & max_squared_distance = [](
    const AABB<DerivedV,DIM> * A,
    const AABB<Derivedother_V,DIM> * B)->Scalar
  {
    AlignedBox<Scalar,DIM> combo = A->m_box;
    combo.extend(B->m_box);
    return combo.diagonal().squaredNorm();
  };

  //// other base-case
  //if(other->is_leaf())
  //{
  //  double sqr_d = sqrD(other->m_primitive);
  //  int i = I(other->m_primitive);
  //  RowVectorDIMS c = C.row(m_primitive);
  //  RowVectorDIMS p = other_V.row(m_primitive);
  //  leaf_squared_distance(V,Ele,p,sqr_d,i,c);
  //  sqrD(other->m_primitive) = sqr_d;
  //  I(other->m_primitive) = i;
  //  C.row(m_primitive) = c;
  //  return;
  //}
  std::vector<const AABB<DerivedV,DIM> * > this_list;
  if(this->is_leaf())
  {
    this_list.push_back(this);
  }else
  {
    assert(this->m_left);
    this_list.push_back(this->m_left);
    assert(this->m_right);
    this_list.push_back(this->m_right);
  }
  std::vector<AABB<Derivedother_V,DIM> *> other_list;
  if(other->is_leaf())
  {
    other_list.push_back(other);
  }else
  {
    assert(other->m_left);
    other_list.push_back(other->m_left);
    assert(other->m_right);
    other_list.push_back(other->m_right);
  }

  //const std::function<Scalar(
  //  const AABB<Derivedother_V,DIM> * other)
  //    > max_sqr_d = [&sqrD,&max_sqr_d](const AABB<Derivedother_V,DIM> * other)->Scalar
  //  {
  //    if(other->is_leaf())
  //    {
  //      return sqrD(other->m_primitive);
  //    }else
  //    {
  //      return std::max(max_sqr_d(other->m_left),max_sqr_d(other->m_right));
  //    }
  //  };

  //// Potentially recurse on all pairs, if minimum distance is less than running
  //// bound
  //Eigen::Matrix<Scalar,Eigen::Dynamic,1> other_max_sqr_d =
  //  Eigen::Matrix<Scalar,Eigen::Dynamic,1>::Constant(other_list.size(),1,min_sqr_d);
  for(size_t child = 0;child<other_list.size();child++)
  {
    auto other_tree = other_list[child];

    Eigen::Matrix<Scalar,Eigen::Dynamic,1> this_max_sqr_d(this_list.size(),1);
    for(size_t t = 0;t<this_list.size();t++)
    {
      const auto this_tree = this_list[t];
      this_max_sqr_d(t) = max_squared_distance(this_tree,other_tree);
    }
    if(this_list.size() ==2 &&
      ( this_max_sqr_d(0) > this_max_sqr_d(1))
      )
    {
      std::swap(this_list[0],this_list[1]);
      //std::swap(this_max_sqr_d(0),this_max_sqr_d(1));
    }
    const Scalar sqr_d = this_max_sqr_d.minCoeff();


    for(size_t t = 0;t<this_list.size();t++)
    {
      const auto this_tree = this_list[t];

      //const auto mm = max_sqr_d(other_tree);
      //const Scalar mc = other_max_sqr_d(child);
      //assert(mc == mm);
      // Only look left/right in this_list if can possible decrease somebody's
      // distance in this_tree.
      const Scalar min_this_other = min_squared_distance(this_tree,other_tree); 
      if(
          min_this_other < sqr_d && 
          min_this_other < other_tree->m_max_sqr_d)
      {
        //cout<<"before: "<<other_max_sqr_d(child)<<endl;
        //other_max_sqr_d(child) = std::min(
        //  other_max_sqr_d(child),
        //  this_tree->squared_distance_helper(
        //    V,Ele,other_tree,other_V,other_Ele,other_max_sqr_d(child),sqrD,I,C));
        //cout<<"after: "<<other_max_sqr_d(child)<<endl;
          this_tree->squared_distance_helper(
            V,Ele,other_tree,other_V,other_Ele,0,sqrD,I,C);
      }
    }
  }
  //const Scalar ret = other_max_sqr_d.maxCoeff();
  //const auto mm = max_sqr_d(other);
  //assert(mm == ret);
  //cout<<"non-leaf: "<<ret<<endl;
  //return ret;
  if(!other->is_leaf())
  {
    other->m_max_sqr_d = std::max(other->m_left->m_max_sqr_d,other->m_right->m_max_sqr_d);
  }
  return 0;
#endif
}

template <typename DerivedV, int DIM>
inline void igl::AABB<DerivedV,DIM>::leaf_squared_distance(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::MatrixXi & Ele, 
  const RowVectorDIMS & p,
  Scalar & sqr_d,
  int & i,
  RowVectorDIMS & c) const
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;

  // Simplex size
  const size_t ss = Ele.cols();
  // Only one element per node
  // plane unit normal
  bool inside_triangle = false;
  Scalar d_j = std::numeric_limits<Scalar>::infinity();
  RowVectorDIMS pp;
  // Only consider triangles, and non-degenerate triangles at that
  if(ss == 3 && 
      Ele(m_primitive,0) != Ele(m_primitive,1) && 
      Ele(m_primitive,1) != Ele(m_primitive,2) && 
      Ele(m_primitive,2) != Ele(m_primitive,0))
  {
    const RowVectorDIMS v10 = (V.row(Ele(m_primitive,1))- V.row(Ele(m_primitive,0)));
    const RowVectorDIMS v20 = (V.row(Ele(m_primitive,2))- V.row(Ele(m_primitive,0)));
    const RowVectorDIMS n = v10.cross(v20);
    Scalar n_norm = n.norm();
    if(n_norm > 0)
    {
      const RowVectorDIMS un = n/n.norm();
      // vector to plane
      const RowVectorDIMS bc = 
        1./3.*
        ( V.row(Ele(m_primitive,0))+
          V.row(Ele(m_primitive,1))+
          V.row(Ele(m_primitive,2)));
      const auto & v = p-bc;
      // projected point on plane
      d_j = v.dot(un);
      pp = p - d_j*un;
      // determine if pp is inside triangle
      Eigen::Matrix<Scalar,1,3> b;
      barycentric_coordinates(
            pp,
            V.row(Ele(m_primitive,0)),
            V.row(Ele(m_primitive,1)),
            V.row(Ele(m_primitive,2)),
            b);
      inside_triangle = fabs(fabs(b(0)) + fabs(b(1)) + fabs(b(2)) - 1.) <= 1e-10;
    }
  }
  const auto & point_point_squared_distance = [&](const RowVectorDIMS & s)
  {
    const Scalar sqr_d_s = (p-s).squaredNorm();
    set_min(p,sqr_d_s,m_primitive,s,sqr_d,i,c);
  };
  if(inside_triangle)
  {
    // point-triangle squared distance
    const Scalar sqr_d_j = d_j*d_j;
    //cout<<"point-triangle..."<<endl;
    set_min(p,sqr_d_j,m_primitive,pp,sqr_d,i,c);
  }else
  {
    if(ss >= 2)
    {
      // point-segment distance
      // number of edges
      size_t ne = ss==3?3:1;
      for(size_t x = 0;x<ne;x++)
      {
        const size_t e1 = Ele(m_primitive,(x+1)%ss);
        const size_t e2 = Ele(m_primitive,(x+2)%ss);
        const RowVectorDIMS & s = V.row(e1);
        const RowVectorDIMS & d = V.row(e2);
        // Degenerate edge
        if(e1 == e2 || (s-d).squaredNorm()==0)
        {
          // only consider once
          if(e1 < e2)
          {
            point_point_squared_distance(s);
          }
          continue;
        }
        Matrix<Scalar,1,1> sqr_d_j_x(1,1);
        Matrix<Scalar,1,1> t(1,1);
        project_to_line_segment(p,s,d,t,sqr_d_j_x);
        const RowVectorDIMS q = s+t(0)*(d-s);
        set_min(p,sqr_d_j_x(0),m_primitive,q,sqr_d,i,c);
      }
    }else
    {
      // So then Ele is just a list of points...
      assert(ss == 1);
      const RowVectorDIMS & s = V.row(Ele(m_primitive,0));
      point_point_squared_distance(s);
    }
  }
}


template <typename DerivedV, int DIM>
inline void igl::AABB<DerivedV,DIM>::set_min(
  const RowVectorDIMS & p,
  const Scalar sqr_d_candidate,
  const int i_candidate,
  const RowVectorDIMS & c_candidate,
  Scalar & sqr_d,
  int & i,
  RowVectorDIMS & c) const
{
#ifndef NDEBUG
  //std::cout<<matlab_format(c_candidate,"c_candidate")<<std::endl;
  const Scalar pc_norm = (p-c_candidate).squaredNorm();
  const Scalar diff = fabs(sqr_d_candidate - pc_norm);
  assert(diff<=1e-10 && "distance should match norm of difference");
#endif
  if(sqr_d_candidate < sqr_d)
  {
    i = i_candidate;
    c = c_candidate;
    sqr_d = sqr_d_candidate;
  }
}


template <typename DerivedV, int DIM>
inline void
igl::AABB<DerivedV,DIM>::barycentric_coordinates(
  const RowVectorDIMS & p, 
  const RowVectorDIMS & a, 
  const RowVectorDIMS & b, 
  const RowVectorDIMS & c,
  Eigen::Matrix<Scalar,1,3> & bary)
{
  // http://gamedev.stackexchange.com/a/23745
  const RowVectorDIMS v0 = b - a;
  const RowVectorDIMS v1 = c - a;
  const RowVectorDIMS v2 = p - a;
  Scalar d00 = v0.dot(v0);
  Scalar d01 = v0.dot(v1);
  Scalar d11 = v1.dot(v1);
  Scalar d20 = v2.dot(v0);
  Scalar d21 = v2.dot(v1);
  Scalar denom = d00 * d11 - d01 * d01;
  bary(1) = (d11 * d20 - d01 * d21) / denom;
  bary(2) = (d00 * d21 - d01 * d20) / denom;
  bary(0) = 1.0f - bary(1) - bary(2);
}

#endif
