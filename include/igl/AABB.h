// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_AABB_H
#define IGL_AABB_H

#include "Hit.h"
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
namespace igl
{
  /// Implementation of semi-general purpose axis-aligned bounding box hierarchy.
  /// The mesh (V,Ele) is stored and managed by the caller and each routine here
  /// simply takes it as references (it better not change between calls).
  ///
  /// It's a little annoying that the Dimension is a template parameter and not
  /// picked up at run time from V. This leads to duplicated code for 2d/3d (up to
  /// dim).
  ///
  /// @tparam DerivedV  Matrix type of vertex positions (e.g., `Eigen::MatrixXd`)
  /// @tparam DIM Dimension of mesh vertex positions (2 or 3)
  template <typename DerivedV, int DIM>
    class AABB 
    {
public:
      /// Scalar type of vertex positions (e.g., `double`)
      typedef typename DerivedV::Scalar Scalar;
      /// Fixed-size (`DIM`) RowVector type using `Scalar`
      typedef Eigen::Matrix<Scalar,1,DIM> RowVectorDIMS;
      /// Fixed-size (`DIM`) (Column)Vector type using `Scalar`
      typedef Eigen::Matrix<Scalar,DIM,1> VectorDIMS;
      /// Fixed-width (`DIM`) Matrix type using `Scalar`
      typedef Eigen::Matrix<Scalar,Eigen::Dynamic,DIM> MatrixXDIMS;
      /// Pointer to "left" child node (`nullptr` if leaf)
      // Shared pointers are slower...
      AABB * m_left; 
      /// Pointer to "right" child node (`nullptr` if leaf)
      AABB * m_right;
      /// Axis-Aligned Bounding Box containing this node
      Eigen::AlignedBox<Scalar,DIM> m_box;
      /// Index of single primitive in this node if full leaf, otherwise -1 for non-leaf
      int m_primitive;
      //Scalar m_low_sqr_d;
      //int m_depth;
      /// @private
      AABB():
        m_left(NULL), m_right(NULL),
        m_box(), m_primitive(-1)
        //m_low_sqr_d(std::numeric_limits<double>::infinity()),
        //m_depth(0)
    {
      static_assert(DerivedV::ColsAtCompileTime == DIM || DerivedV::ColsAtCompileTime == Eigen::Dynamic,"DerivedV::ColsAtCompileTime == DIM || DerivedV::ColsAtCompileTime == Eigen::Dynamic");
    }
      /// @private
      // http://stackoverflow.com/a/3279550/148668
      AABB(const AABB& other):
        m_left(other.m_left ? new AABB(*other.m_left) : NULL),
        m_right(other.m_right ? new AABB(*other.m_right) : NULL),
        m_box(other.m_box),
        m_primitive(other.m_primitive)
        //m_low_sqr_d(other.m_low_sqr_d),
        //m_depth(std::max(
        //   m_left ? m_left->m_depth + 1 : 0,
        //   m_right ? m_right->m_depth + 1 : 0))
        {
        }
      /// @private
      // copy-swap idiom
      friend void swap(AABB& first, AABB& second)
      {
        // Enable ADL
        using std::swap;
        swap(first.m_left,second.m_left);
        swap(first.m_right,second.m_right);
        swap(first.m_box,second.m_box);
        swap(first.m_primitive,second.m_primitive);
        //swap(first.m_low_sqr_d,second.m_low_sqr_d);
        //swap(first.m_depth,second.m_depth);
      }
      /// @private
      // Pass-by-value (aka copy)
      AABB& operator=(AABB other)
      {
        swap(*this,other);
        return *this;
      }
      /// @private
      AABB(AABB&& other):
        // initialize via default constructor
        AABB() 
      {
        swap(*this,other);
      }
      /// @private
      // Seems like there should have been an elegant solution to this using
      // the copy-swap idiom above:
      IGL_INLINE void deinit()
      {
        m_primitive = -1;
        m_box = Eigen::AlignedBox<Scalar,DIM>();
        delete m_left;
        m_left = NULL;
        delete m_right;
        m_right = NULL;
      }
      /// @private
      ~AABB()
      {
        deinit();
      }
      /// Build an Axis-Aligned Bounding Box tree for a given mesh and given
      /// serialization of a previous AABB tree.
      ///
      /// @param[in] V  #V by dim list of mesh vertex positions. 
      /// @param[in] Ele  #Ele by dim+1 list of mesh indices into #V. 
      /// @param[in] bb_mins  max_tree by dim list of bounding box min corner positions
      /// @param[in] bb_maxs  max_tree by dim list of bounding box max corner positions
      /// @param[in] elements  max_tree list of element or (not leaf id) indices into Ele
      /// @param[in] i  recursive call index {0}
      template <
        typename DerivedEle, 
        typename Derivedbb_mins, 
        typename Derivedbb_maxs,
        typename Derivedelements>
        IGL_INLINE void init(
            const Eigen::MatrixBase<DerivedV> & V,
            const Eigen::MatrixBase<DerivedEle> & Ele, 
            const Eigen::MatrixBase<Derivedbb_mins> & bb_mins,
            const Eigen::MatrixBase<Derivedbb_maxs> & bb_maxs,
            const Eigen::MatrixBase<Derivedelements> & elements,
            const int i = 0);
      /// Build an Axis-Aligned Bounding Box tree for a given mesh and given
      /// serialization of a previous AABB tree.
      ///
      /// @param[in] V  #V by dim list of mesh vertex positions. 
      /// @param[in] Ele  #Ele by dim+1 list of mesh indices into #V. 
      template <typename DerivedEle>
      IGL_INLINE void init(
          const Eigen::MatrixBase<DerivedV> & V,
          const Eigen::MatrixBase<DerivedEle> & Ele);
      /// Build an Axis-Aligned Bounding Box tree for a given mesh.
      ///
      /// @param[in] V  #V by dim list of mesh vertex positions. 
      /// @param[in] Ele  #Ele by dim+1 list of mesh indices into #V. 
      /// @param[in] SI  #Ele by dim list revealing for each coordinate where Ele's
      ///              barycenters would be sorted: SI(e,d) = i --> the dth coordinate of
      ///              the barycenter of the eth element would be placed at position i in a
      ///              sorted list.
      /// @param[in] I  #I list of indices into Ele of elements to include (for recursive
      ///     calls)
      /// 
      template <typename DerivedEle, typename DerivedSI, typename DerivedI>
      IGL_INLINE void init(
          const Eigen::MatrixBase<DerivedV> & V,
          const Eigen::MatrixBase<DerivedEle> & Ele, 
          const Eigen::MatrixBase<DerivedSI> & SI,
          const Eigen::MatrixBase<DerivedI>& I);
      /// Return whether at leaf node
      IGL_INLINE bool is_leaf() const;
      /// Find the indices of elements containing given point: this makes sense
      /// when Ele is a co-dimension 0 simplex (tets in 3D, triangles in 2D).
      ///
      /// @param[in]  V  #V by dim list of mesh vertex positions. **Should be same as used to
      ///     construct mesh.**
      /// @param[in]  Ele  #Ele by dim+1 list of mesh indices into #V. **Should be same as used to
      ///     construct mesh.**
      /// @param[in]  q  dim row-vector query position
      /// @param[in]  first  whether to only return first element containing q
      /// @return  list of indices of elements containing q
      template <typename DerivedEle, typename Derivedq>
      IGL_INLINE std::vector<int> find(
          const Eigen::MatrixBase<DerivedV> & V,
          const Eigen::MatrixBase<DerivedEle> & Ele, 
          const Eigen::MatrixBase<Derivedq> & q,
          const bool first=false) const;

      /// Number of nodes contained in subtree
      ///
      /// @return Number of elements m then total tree size should be 2*h where h is
      /// the deepest depth 2^ceil(log(#Ele*2-1))
      IGL_INLINE int subtree_size() const;

      /// Serialize this class into 3 arrays (so we can pass it pack to matlab)
      ///
      /// @param[out]  bb_mins  max_tree by dim list of bounding box min corner positions
      /// @param[out]  bb_maxs  max_tree by dim list of bounding box max corner positions
      /// @param[out]  elements  max_tree list of element or (not leaf id) indices into Ele
      /// @param[in]  i  recursive call index into these arrays {0}
      template <
        typename Derivedbb_mins, 
        typename Derivedbb_maxs,
        typename Derivedelements>
        IGL_INLINE void serialize(
            Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
            Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
            Eigen::PlainObjectBase<Derivedelements> & elements,
            const int i = 0) const;
      /// Compute squared distance to a query point
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  p  dim-long query point 
      /// @param[out]  i  facet index corresponding to smallest distances
      /// @param[out]  c  closest point
      /// @return squared distance
      ///
      /// \pre Currently assumes Elements are triangles regardless of
      /// dimension.
      template <typename DerivedEle>
      IGL_INLINE Scalar squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & p,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;
      /// Compute squared distance to a query point if within `low_sqr_d` and
      /// `up_sqr_d`.
      ///
      /// @param[in] V  #V by dim list of vertex positions
      /// @param[in] Ele  #Ele by dim list of simplex indices
      /// @param[in] p  dim-long query point 
      /// @param[in] low_sqr_d  lower bound on squared distance, specified maximum squared
      ///              distance 
      /// @param[in] up_sqr_d  current upper bounded on squared distance, current minimum
      ///              squared distance (only consider distances less than this), see
      ///              output.
      /// @param[out]  i  facet index corresponding to smallest distances
      /// @param[out]  c  closest point
      /// @return squared distance
      ///
      /// \pre currently assumes Elements are triangles regardless of
      /// dimension.
      template <typename DerivedEle>
      IGL_INLINE Scalar squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & p,
        const Scalar low_sqr_d,
        const Scalar up_sqr_d,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;
      /// Compute squared distance to a query point (default `low_sqr_d`)
      ///
      /// @param[in] V  #V by dim list of vertex positions
      /// @param[in] Ele  #Ele by dim list of simplex indices
      /// @param[in] p  dim-long query point 
      /// @param[in] up_sqr_d  current upper bounded on squared distance, current minimum
      ///              squared distance (only consider distances less than this), see
      ///              output.
      /// @param[out]  i  facet index corresponding to smallest distances
      /// @param[out]  c  closest point
      /// @return squared distance
      ///
      template <typename DerivedEle>
      IGL_INLINE Scalar squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & p,
        const Scalar up_sqr_d,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;
      /// Intersect a ray with the mesh return all hits
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  origin  dim-long ray origin
      /// @param[in]  dir  dim-long ray direction
      /// @param[out]  hits  list of hits
      /// @return  true if any hits
      template <typename DerivedEle>
      IGL_INLINE bool intersect_ray(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & origin,
        const RowVectorDIMS & dir,
        std::vector<igl::Hit> & hits) const;
      /// Intersect a ray with the mesh return first hit
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  origin  dim-long ray origin
      /// @param[in]  dir  dim-long ray direction
      /// @param[out]  hit  first hit
      /// @return  true if any hit
      template <typename DerivedEle>
      IGL_INLINE bool intersect_ray(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & origin,
        const RowVectorDIMS & dir,
        igl::Hit & hit) const;
      /// Intersect a ray with the mesh return first hit farther than `min_t`
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  origin  dim-long ray origin
      /// @param[in]  dir  dim-long ray direction
      /// @param[in]  min_t  minimum t value to consider
      /// @param[out]  hit  first hit
      /// @return  true if any hit
      template <typename DerivedEle>
      IGL_INLINE bool intersect_ray(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & origin,
        const RowVectorDIMS & dir,
        const Scalar min_t,
        igl::Hit & hit) const;

public:
      /// Compute the squared distance from all query points in P to the
      /// _closest_ points on the primitives stored in the AABB hierarchy for
      /// the mesh (V,Ele).
      ///
      /// @param[in] V  #V by dim list of vertex positions
      /// @param[in] Ele  #Ele by dim list of simplex indices
      /// @param[in] P  #P by dim list of query points
      /// @param[out] sqrD  #P list of squared distances
      /// @param[out] I  #P list of indices into Ele of closest primitives
      /// @param[out] C  #P by dim list of closest points
      template <
        typename DerivedEle,
        typename DerivedP, 
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      IGL_INLINE void squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const Eigen::MatrixBase<DerivedP> & P,
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;

      /// Compute the squared distance from all query points in P already stored
      /// in its own AABB hierarchy to the _closest_ points on the primitives
      /// stored in the AABB hierarchy for the mesh (V,Ele).
      ///
      /// @param[in] V  #V by dim list of vertex positions
      /// @param[in] Ele  #Ele by dim list of simplex indices
      /// @param[in] other  AABB hierarchy of another set of primitives (must be points)
      /// @param[in] other_V  #other_V by dim list of query points
      /// @param[in] other_Ele  #other_Ele by ss list of simplex indices into other_V
      ///     (must be simple list of points: ss == 1)
      /// @param[out] sqrD  #P list of squared distances
      /// @param[out] I  #P list of indices into Ele of closest primitives
      /// @param[out] C  #P by dim list of closest points
      template < 
        typename DerivedEle,
        typename Derivedother_V,
        typename Derivedother_Ele,
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      IGL_INLINE void squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const AABB<Derivedother_V,DIM> & other,
        const Eigen::MatrixBase<Derivedother_V> & other_V,
        const Eigen::MatrixBase<Derivedother_Ele> & other_Ele, 
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;
private:
      template < 
        typename DerivedEle,
        typename Derivedother_V,
        typename Derivedother_Ele,
        typename DerivedsqrD, 
        typename DerivedI, 
        typename DerivedC>
      IGL_INLINE Scalar squared_distance_helper(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const AABB<Derivedother_V,DIM> * other,
        const Eigen::MatrixBase<Derivedother_V> & other_V,
        const Eigen::MatrixBase<Derivedother_Ele>& other_Ele, 
        const Scalar up_sqr_d,
        Eigen::PlainObjectBase<DerivedsqrD> & sqrD,
        Eigen::PlainObjectBase<DerivedI> & I,
        Eigen::PlainObjectBase<DerivedC> & C) const;
      // Compute the squared distance to the primitive in this node: assumes
      // that this is indeed a leaf node.
      //
      // Inputs:
      //   V  #V by dim list of vertex positions
      //   Ele  #Ele by dim list of simplex indices
      //   p  dim-long query point
      //   sqr_d  current minimum distance for this query, see output
      //   i  current index into Ele of closest point, see output
      //   c  dim-long current closest point, see output
      // Outputs:
      //   sqr_d   minimum of initial value and squared distance to this
      //     primitive
      //   i  possibly updated index into Ele of closest point
      //   c  dim-long possibly updated closest point
      template <typename DerivedEle>
      IGL_INLINE void leaf_squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & p,
        const Scalar low_sqr_d,
        Scalar & sqr_d,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;
      // Default low_sqr_d
      template <typename DerivedEle>
      IGL_INLINE void leaf_squared_distance(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele, 
        const RowVectorDIMS & p,
        Scalar & sqr_d,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;
      // If new distance (sqr_d_candidate) is less than current distance
      // (sqr_d), then update this distance and its associated values
      // _in-place_:
      //
      // Inputs:
      //   p  dim-long query point (only used in DEBUG mode)
      //   sqr_d  candidate minimum distance for this query, see output
      //   i  candidate index into Ele of closest point, see output
      //   c  dim-long candidate closest point, see output
      //   sqr_d  current minimum distance for this query, see output
      //   i  current index into Ele of closest point, see output
      //   c  dim-long current closest point, see output
      // Outputs:
      //   sqr_d   minimum of initial value and squared distance to this
      //     primitive
      //   i  possibly updated index into Ele of closest point
      //   c  dim-long possibly updated closest point
      IGL_INLINE void set_min(
        const RowVectorDIMS & p,
        const Scalar sqr_d_candidate,
        const int i_candidate,
        const RowVectorDIMS & c_candidate,
        Scalar & sqr_d,
        int & i,
        Eigen::PlainObjectBase<RowVectorDIMS> & c) const;

      /// Intersect a ray with the mesh return all hits
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  origin  dim-long ray origin
      /// @param[in]  dir  dim-long ray direction
      /// @param[out]  hits  list of hits
      /// @return  true if any hits
      template <typename DerivedEle>
      IGL_INLINE bool intersect_ray_opt(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele,
        const RowVectorDIMS & origin,
        const RowVectorDIMS & dir,
        const RowVectorDIMS & inv_dir,
        const RowVectorDIMS & inv_dir_pad,
        std::vector<igl::Hit> & hits) const;
      /// Intersect a ray with the mesh return first hit farther than `min_t`
      ///
      /// @param[in]  V  #V by dim list of vertex positions
      /// @param[in]  Ele  #Ele by dim list of simplex indices
      /// @param[in]  origin  dim-long ray origin
      /// @param[in]  dir  dim-long ray direction
      /// @param[in]  min_t  minimum t value to consider
      /// @param[out]  hit  first hit
      /// @return  true if any hit
      template <typename DerivedEle>
      IGL_INLINE bool intersect_ray_opt(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedEle> & Ele,
        const RowVectorDIMS & origin,
        const RowVectorDIMS & dir,
        const RowVectorDIMS & inv_dir,
        const RowVectorDIMS & inv_dir_pad,
        const Scalar min_t,
        igl::Hit & hit) const;
public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
}


#ifndef IGL_STATIC_LIBRARY
#  include "AABB.cpp"
#endif

#endif
