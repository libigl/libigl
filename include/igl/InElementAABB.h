#ifndef IGL_IN_ELEMENT_AABB_H
#define IGL_IN_ELEMENT_AABB_H

#include <Eigen/Core>
#include <memory>
#include <vector>
namespace igl
{
  class InElementAABB
  {
    public:
      std::shared_ptr<InElementAABB> m_left, m_right;
      Eigen::RowVectorXd m_bb_min,m_bb_max;
      // -1 non-leaf
      int m_element;
      InElementAABB():
        m_left(NULL), m_right(NULL),
        m_bb_min(), m_bb_max(), m_element(-1)
      {}
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
      inline void init(
          const Eigen::MatrixXd & V,
          const Eigen::MatrixXi & Ele, 
          const Eigen::MatrixXd & bb_mins,
          const Eigen::MatrixXd & bb_maxs,
          const Eigen::VectorXi & elements,
          const int i = 0);
      // Wrapper for root with empty serialization
      inline void init(
        const Eigen::MatrixXd & V,
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
          const Eigen::MatrixXd & V,
          const Eigen::MatrixXi & Ele, 
          const Eigen::MatrixXi & SI,
          const Eigen::VectorXi & I);
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
      inline std::vector<int> find(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & Ele, 
        const Eigen::RowVectorXd & q, 
        const bool first=false) const;
  
      // If number of elements m then total tree size should be 2*h where h is
      // the deepest depth 2^ceil(log(#Ele*2-1))
      inline int subtree_size();
  
      // Serialize this class into 3 arrays (so we can pass it pack to matlab)
      //
      // Outputs:
      //   bb_mins  max_tree by dim list of bounding box min corner positions
      //   bb_maxs  max_tree by dim list of bounding box max corner positions
      //   elements  max_tree list of element or (not leaf id) indices into Ele
      //   i  recursive call index into these arrays {0}
      inline void serialize(
        Eigen::MatrixXd & bb_mins,
        Eigen::MatrixXd & bb_maxs,
        Eigen::VectorXi & elements,
        const int i = 0);
  };
}

// Implementation
#include <igl/volume.h>
#include <igl/colon.h>
#include <igl/doublearea.h>
#include <igl/matlab_format.h>
#include <igl/colon.h>
#include <igl/sort.h>
#include <igl/barycenter.h>
#include <iostream>

inline void igl::InElementAABB::init(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::MatrixXd & bb_mins,
    const Eigen::MatrixXd & bb_maxs,
    const Eigen::VectorXi & elements,
    const int i)
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  if(bb_mins.size() > 0)
  {
    assert(bb_mins.rows() == bb_maxs.rows() && "Serial tree arrays must match");
    assert(bb_mins.cols() == V.cols() && "Serial tree array dim must match V");
    assert(bb_mins.cols() == bb_maxs.cols() && "Serial tree arrays must match");
    assert(bb_mins.rows() == elements.rows() &&
      "Serial tree arrays must match");
    // construct from serialization
    m_bb_min = bb_mins.row(i);
    m_bb_max = bb_maxs.row(i);
    m_element = elements(i);
    // Not leaf then recurse
    if(m_element == -1)
    {
      m_left = make_shared<InElementAABB>();
      m_left->init( V,Ele,bb_mins,bb_maxs,elements,2*i+1);
      m_right = make_shared<InElementAABB>();
      m_right->init( V,Ele,bb_mins,bb_maxs,elements,2*i+2);
    }
  }else
  {
    VectorXi allI = colon<int>(0,Ele.rows()-1);
    MatrixXd BC;
    barycenter(V,Ele,BC);
    MatrixXi SI(BC.rows(),BC.cols());
    {
      MatrixXd _;
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

inline void igl::InElementAABB::init(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele)
{
  using namespace Eigen;
  return init(V,Ele,MatrixXd(),MatrixXd(),VectorXi(),0);
}

inline void igl::InElementAABB::init(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::MatrixXi & SI,
    const Eigen::VectorXi & I)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  const int dim = V.cols();
  const double inf = numeric_limits<double>::infinity();
  m_bb_min.setConstant(1,dim, inf);
  m_bb_max.setConstant(1,dim,-inf);
  // Compute bounding box
  for(int i = 0;i<I.rows();i++)
  {
    for(int c = 0;c<Ele.cols();c++)
    {
      for(int d = 0;d<dim;d++)
      {
        m_bb_min(d) = min(m_bb_min(d),V(Ele(I(i),c),d));
        m_bb_max(d) = max(m_bb_max(d),V(Ele(I(i),c),d));
      }
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
        m_element = I(0);
        break;
      }
    default:
      {
        // Compute longest direction
        int max_d = -1;
        double max_len = -inf;
        for(int d = 0;d<dim;d++)
        {
          const auto diff = (m_bb_max[d] - m_bb_min[d]);
          if( diff > max_len )
          {
            max_len = diff;
            max_d = d;
          }
        }
        // Can't use median on BC directly because many may have same value,
        // but can use median on sorted BC indices
        VectorXi SIdI(I.rows());
        for(int i = 0;i<I.rows();i++)
        {
          SIdI(i) = SI(I(i),max_d);
        }
        // Since later I use <= I think I don't need to worry about odd/even
        // Pass by copy to avoid changing input
        const auto median = [](VectorXi A)->double
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
        const double med = median(SIdI);
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
        if(LI.rows()>0)
        {
          m_left = make_shared<InElementAABB>();
          m_left->init(V,Ele,SI,LI);
        }
        if(RI.rows()>0)
        {
          m_right = make_shared<InElementAABB>();
          m_right->init(V,Ele,SI,RI);
        }
      }
  }
}

inline std::vector<int> igl::InElementAABB::find(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele, 
    const Eigen::RowVectorXd & q, 
    const bool first) const
{
  using namespace std;
  using namespace igl;
  using namespace Eigen;
  bool inside = true;
  const int dim = m_bb_max.size();
  assert(q.size() == m_bb_max.size());
  const double epsilon = 1e-14;
  // Check if outside bounding box
  for(int d = 0;d<q.size()&&inside;d++)
  {
    inside &= (q(d)-m_bb_min(d))>=epsilon;
    inside &= (m_bb_max(d)-q(d))>=epsilon;
  }
  cout<<"searching..."<<endl;
  if(!inside)
  {
    cout<<"not in bb"<<endl;
    return std::vector<int>();
  }
  if(m_element != -1)
  {
    // Initialize to some value > -epsilon
    double a1=1,a2=1,a3=1,a4=1;
    switch(dim)
    {
      case 3:
        {
          // Barycentric coordinates
          const RowVector3d V1 = V.row(Ele(m_element,0));
          const RowVector3d V2 = V.row(Ele(m_element,1));
          const RowVector3d V3 = V.row(Ele(m_element,2));
          const RowVector3d V4 = V.row(Ele(m_element,3));
          a1 = volume_single(V2,V4,V3,(RowVector3d)q);
          a2 = volume_single(V1,V3,V4,(RowVector3d)q);
          a3 = volume_single(V1,V4,V2,(RowVector3d)q);
          a4 = volume_single(V1,V2,V3,(RowVector3d)q);
          break;
        }
      case 2:
        {
          // Barycentric coordinates
          const Vector2d V1 = V.row(Ele(m_element,0));
          const Vector2d V2 = V.row(Ele(m_element,1));
          const Vector2d V3 = V.row(Ele(m_element,2));
          double a0 = doublearea_single(V1,V2,V3);
          a1 = doublearea_single(V1,V2,(Vector2d)q);
          a2 = doublearea_single(V2,V3,(Vector2d)q);
          a3 = doublearea_single(V3,V1,(Vector2d)q);
          cout<<
            a0<<" "<<
            a1<<" "<<
            a2<<" "<<
            a3<<" "<<
            endl;
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
      return std::vector<int>(1,m_element);
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

inline int igl::InElementAABB::subtree_size()
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


inline void igl::InElementAABB::serialize(
    Eigen::MatrixXd & bb_mins,
    Eigen::MatrixXd & bb_maxs,
    Eigen::VectorXi & elements,
    const int i)
{
  using namespace std;
  // Calling for root then resize output
  if(i==0)
  {
    const int m = subtree_size();
    //cout<<"m: "<<m<<endl;
    bb_mins.resize(m,m_bb_min.size());
    bb_maxs.resize(m,m_bb_max.size());
    elements.resize(m,1);
  }
  //cout<<i<<" ";
  bb_mins.row(i) = m_bb_min;
  bb_maxs.row(i) = m_bb_max;
  elements(i) = m_element;
  if(m_left != NULL)
  {
    m_left->serialize(bb_mins,bb_maxs,elements,2*i+1);
  }
  if(m_right != NULL)
  {
    m_right->serialize(bb_mins,bb_maxs,elements,2*i+2);
  }
}
#endif
