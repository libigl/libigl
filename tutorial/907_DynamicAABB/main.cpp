#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/unique_simplices.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/barycenter.h>
#include <igl/writePLY.h>
#include <igl/STR.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
#include <igl/writeDMAT.h>
#include <igl/decimate.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/remove_unreferenced.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/is_self_intersecting.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/edge_flaps.h>
#include <igl/doublearea.h>
#include <igl/decimate_callback_types.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/tri_tri_intersect.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/PI.h>
#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/matlab_format.h>
#include <limits>
#include <igl/randperm.h>
#include <igl/avg_edge_length.h>
#include <igl/find.h>
#include <deque>
extern "C"
{
#include <igl/raytri.c>
}





const int MAX_RUNS = 10;
using AABB = igl::AABB<Eigen::MatrixXd,3>;

template <typename DerivedV, int DIM>
void validate(const igl::AABB<DerivedV,DIM> * tree, int depth = 0)
{
  if(tree->is_leaf())
  {
    assert(tree->m_primitive >= 0 || tree->is_root());
  }
  if(tree->m_left)
  {
    assert(tree->m_box.contains(tree->m_left->m_box));
    assert(tree->m_left->m_parent == tree);
    validate(tree->m_left,depth+1);
  }
  if(tree->m_right)
  {
    assert(tree->m_box.contains(tree->m_right->m_box));
    assert(tree->m_right->m_parent == tree);
    validate(tree->m_right,depth+1);
  }
}

template <typename DerivedV, int DIM>
void validate(
    const igl::AABB<DerivedV,DIM> * root,
    const std::vector<igl::AABB<DerivedV,DIM> > & leaves)
{
  // Check that all leaves are in the tree
  for(int i = 0;i<leaves.size();i++)
  {
    const auto * leaf = &leaves[i];
    assert(leaf->m_primitive == i);
    assert(leaf->root() == root);
  }
  return validate(root);
}

template <typename DerivedV, int DIM>
void print(const igl::AABB<DerivedV,DIM> * tree, int depth = 0)
{
  const auto indent = std::string(depth*2,' ');
  printf("%s%p",indent.c_str(),tree);
  if(tree->is_leaf())
  {
    printf(" [%d]",tree->m_primitive);
  }
  printf("\n");
  if(tree->m_left)
  {
    assert(tree->m_box.contains(tree->m_left->m_box));
    assert(tree->m_left->m_parent == tree);
    print(tree->m_left,depth+1);
  }
  if(tree->m_right)
  {
    assert(tree->m_box.contains(tree->m_right->m_box));
    assert(tree->m_right->m_parent == tree);
    print(tree->m_right,depth+1);
  }
}


template <typename DerivedV, int DIM>
int size(const igl::AABB<DerivedV,DIM> * tree)
{
  return 1 + (tree->m_left?size(tree->m_left):0) + (tree->m_right?size(tree->m_right):0);
}
template <typename DerivedV, int DIM>
int height(const igl::AABB<DerivedV,DIM> * tree)
{
  return 1 + std::max((tree->m_left?height(tree->m_left):0),(tree->m_right?height(tree->m_right):0));
}

void box_faces(
  const Eigen::AlignedBox<double,3> & box,
  const double shrink,
  Eigen::MatrixXd & P,
  Eigen::MatrixXi & Q)
{
  auto min_corner = box.min();
  auto max_corner = box.max();
  // shrink by 3%
  min_corner = min_corner + shrink*(max_corner-min_corner);
  max_corner = max_corner - shrink*(max_corner-min_corner);
  P.resize(8,3);
  Q.resize(6,4);
  int p = 0;
  int q = 0;
  Q.row(q++) << p+0,p+1,p+2,p+3;
  Q.row(q++) << p+0,p+1,p+5,p+4;
  Q.row(q++) << p+1,p+2,p+6,p+5;
  Q.row(q++) << p+2,p+3,p+7,p+6;
  Q.row(q++) << p+3,p+0,p+4,p+7;
  Q.row(q++) << p+4,p+5,p+6,p+7;
  P.row(p++) = min_corner;
  P.row(p++) = Eigen::RowVector3d(max_corner[0],min_corner[1],min_corner[2]);
  P.row(p++) = Eigen::RowVector3d(max_corner[0],max_corner[1],min_corner[2]);
  P.row(p++) = Eigen::RowVector3d(min_corner[0],max_corner[1],min_corner[2]);
  P.row(p++) = Eigen::RowVector3d(min_corner[0],min_corner[1],max_corner[2]);
  P.row(p++) = Eigen::RowVector3d(max_corner[0],min_corner[1],max_corner[2]);
  P.row(p++) = max_corner;
  P.row(p++) = Eigen::RowVector3d(min_corner[0],max_corner[1],max_corner[2]);
}

template <typename DerivedV, int DIM>
void box_faces(
  const igl::AABB<DerivedV,DIM> & tree,
  Eigen::MatrixXd & P,
  Eigen::MatrixXi & Q,
  Eigen::VectorXi & D)
{
  static_assert(DIM == 3,"Assumes 3D");
  const int num_nodes = size(&tree);
  P.resize(8*num_nodes,DIM);
  Q.resize(6*num_nodes,4);
  D.resize(6*num_nodes);
  int d = 0;
  int p = 0;
  int q = 0;
  std::vector<std::pair<const igl::AABB<DerivedV,DIM> *,int> > stack;
  stack.push_back({&tree,0});
  while(!stack.empty())
  {
    const auto pair = stack.back();
    const auto * node = pair.first;
    const int depth = pair.second;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    stack.pop_back();
    const auto & box = node->m_box;

    Eigen::MatrixXd Pi;
    Eigen::MatrixXi Qi;
    box_faces(box,0.03,Pi,Qi);
    P.block(p,0,8,DIM) = Pi;
    Q.block(q,0,6,4) = Qi.array()+p;
    p += 8;
    q += 6;
    if(node->m_left)
    {
      stack.push_back({node->m_left,depth+1});
    }
    if(node->m_right)
    {
      stack.push_back({node->m_right,depth+1});
    }
  }
}

void quad_edges(const Eigen::MatrixXi & Q, Eigen::MatrixXi & E)
{
  E.resize(4*Q.rows(),2);
  E <<
    Q.col(0), Q.col(1),
    Q.col(1), Q.col(2),
    Q.col(2), Q.col(3),
    Q.col(3), Q.col(0);
  igl::unique_simplices(Eigen::MatrixXi(E),E);
}

void vis(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const AABB & tree)
{
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
  vr.data().set_face_based(true);
  int depth = 0;
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TQ;
  Eigen::VectorXi TD;
  box_faces(tree,TV,TQ,TD);
  const auto update_edges = [&]()
  {
    Eigen::MatrixXi TQd = TQ(igl::find((TD.array()==depth).eval()),Eigen::all);
    Eigen::MatrixXi TE;
    quad_edges(TQd,TE);
    //Eigen::MatrixXi TF(TQd.rows()*2,3);
    //TF<< 
    //  TQd.col(0),TQd.col(1),TQd.col(2),
    //  TQd.col(0),TQd.col(2),TQd.col(3);
    //vr.append_mesh();
    //vr.data().set_mesh(TV,TF);
    vr.data().set_edges(TV,TE,Eigen::RowVector3d(1,1,1));
    vr.data().line_width = 2;
  };
  update_edges();
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ',':
      case '.':
        depth = std::max(depth+(key==','?-1:1),0);
        update_edges();
        return true;
      default:
        return false;
    }
  };
  vr.launch();
}


/// Determine whether two triangles intersect. We consider the `f`th and `g`th
/// triangles in `F` indexing rows of `V` for 3D positions, but the `c`th corner
/// of the `f`th triangle is replaced by `p`. In matlab, this would be
///
/// ```matlab
/// Tf = V(F(f,:),:);
/// Tf(c,:) = p;
/// ```
///
/// and 
///
/// ```matlab
/// Tg = V(F(g,:),:);
/// ```
///
/// Triangles can share an edge, but only if it's the one opposite the replaced
/// corner. 
///
/// @param[in] V  #V by 3 list of vertex positions
/// @param[in] F  #F by 3 list of triangle indices into rows of V
/// @param[in] E #E by 2 list of unique undirected edge indices into rows of V
/// @param[in] EMAP  #F*3 list of indices into F, mapping each directed edge to
///   unique edge in {1,...,E}
/// @param[in] EF  #E by 2 list of edge indices into F
/// @param[in] EI  #E by 2 list of edge indices into V
/// @param[in] f  index into F of first triangle
/// @param[in] c  index into F of corner of first triangle to replace with `p`
/// @param[in] p  3D position to replace corner of first triangle
/// @param[in] g  index into F of second triangle
/// @returns  true if triangles intersect
///
/// \seealso edge_flaps
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedE,
  typename DerivedEMAP,
  typename DerivedEF,
  typename DerivedEI,
  typename Derivedp>
bool triangle_triangle_intersect(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedEMAP> & EMAP,
  const Eigen::MatrixBase<DerivedEI> & EF,
  const Eigen::MatrixBase<DerivedEF> & EI,
  const int f,
  const int c,
  const Eigen::MatrixBase<Derivedp> & p,
  const int g)
{
  constexpr bool stinker = false;
  if(stinker) { printf("👀\n"); }
  bool found_intersection = false;
  // So edge opposite F(f,c) is the outer edge.
  const int o = EMAP(f + c*F.rows());
  // Do they share the edge opposite c?
  if((EF(o,0) == f && EF(o,1) == g) || (EF(o,1) == f && EF(o,0) == g))
  {
    if(stinker) { printf("⚠️ shares an edge\n"); }
    // Only intersects if the dihedral angle is zero (precondition: no zero
    // area triangles before or after collapse)
#       warning "Consider maintaining a list of (area-vectors)normals?"
    const auto vg10 = (V.row(F(g,1))-V.row(F(g,0))).template head<3>();
    const auto vg20 = (V.row(F(g,2))-V.row(F(g,0))).template head<3>();
    const auto ng = vg10.cross(vg20);
    const int fo = EF(o,0) == f ? EI(o,0) : EI(o,1);
    const auto vf1p = (V.row(F(f,(fo+1)%3))-p).template head<3>();
    const auto vf2p = (V.row(F(f,(fo+2)%3))-p).template head<3>();
    const auto nf = vf1p.cross(vf2p);
    const auto o_vec_un = (V.row(E(o,1))-V.row(E(o,0))).template head<3>();
    const auto o_vec = o_vec_un.stableNormalized();

    const auto dihedral_angle = igl::PI - std::atan2(o_vec.dot(ng.cross(nf)),ng.dot(nf));
    if(dihedral_angle > 1e-8)
    {
      return false;
    }
    // Triangles really really might intersect.
    found_intersection = true;
  }else
  {
    if(stinker) { printf("does not share an edge\n"); }
    // Do they share a vertex?
    int sf,sg;
    bool found_shared_vertex = false;
    for(sf = 0;sf<3;sf++)
    {
      if(sf == c){ continue;}
      for(sg = 0;sg<3;sg++)
      {
        if(F(f,sf) == F(g,sg))
        {
          found_shared_vertex = true;
          break;
        }
      }
      if(found_shared_vertex) { break;} 
    }
    if(found_shared_vertex)
    {
      if(stinker) { printf("⚠️ shares an vertex\n"); }
      // If they share a vertex and intersect, then an opposite edge must
      // stab through the other triangle.

      // intersect_triangle1 needs non-const inputs.
      Eigen::RowVector3d g0 = V.row(F(g,0));
      Eigen::RowVector3d g1 = V.row(F(g,1));
      Eigen::RowVector3d g2 = V.row(F(g,2));
      Eigen::RowVector3d fs;
      if(((sf+1)%3)==c)
      {
        fs = p;
      }else
      {
        fs = V.row(F(f,(sf+1)%3));
      }
      Eigen::RowVector3d fd;
      if( ((sf+2)%3)==c )
      {
        fd = p - fs;
      }else
      {
        fd = V.row(F(f,(sf+2)%3)) - fs;
      }
      double t,u,v;

      if(intersect_triangle1(
            fs.data(),fd.data(),
            g0.data(),g1.data(),g2.data(),
            &t,&u,&v))
      {
        found_intersection = t > 0 && t<1+1e-8;
      }
      if(!found_intersection)
      {
        Eigen::RowVector3d fv[3];
        fv[0] = V.row(F(f,0));
        fv[1] = V.row(F(f,1));
        fv[2] = V.row(F(f,2));
        fv[c] = p;
        Eigen::RowVector3d gs = V.row(F(g,(sg+1)%3));
        Eigen::RowVector3d gd = V.row(F(g,(sg+2)%3)) - gs;
        if(intersect_triangle1(
              gs.data(),gd.data(),
              fv[0].data(),fv[1].data(),fv[2].data(),
              &t,&u,&v))
        {
          found_intersection = t > 0 && t<1+1e-8;
        }
      }
    }else
    {
      bool coplanar;
      Eigen::RowVector3d i1,i2;
      found_intersection = igl::tri_tri_intersection_test_3d(
          V.row(F(g,0)), V.row(F(g,1)), V.row(F(g,2)),
          p,V.row(F(f,(c+1)%3)),V.row(F(f,(c+2)%3)),
          coplanar,
          i1,i2);
      if(stinker) { printf("tri_tri_intersection_test_3d says %s\n",found_intersection?"☠️":"✅"); }
    }
  }
    if(stinker) { printf("%s\n",found_intersection?"☠️":"✅"); }
  return found_intersection;
}

bool collapse_edge_would_create_intersections(
  const int e,
  const Eigen::RowVectorXd & p,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  const igl::AABB<Eigen::MatrixXd,3> & tree)
{
  // Merge two lists of integers
  const auto merge = [&](
    const std::vector<int> & A, const std::vector<int> & B)->
    std::vector<int>
  {
    std::vector<int> C;
    C.reserve( A.size() + B.size() ); // preallocate memory
    C.insert( C.end(), A.begin(), A.end() );
    C.insert( C.end(), B.begin(), B.end() );
    // https://stackoverflow.com/a/1041939/148668
    std::sort( C.begin(), C.end() );
    C.erase( std::unique( C.begin(), C.end() ), C.end() );
    return C;
  };

  std::vector<int> old_one_ring;
  {
    std::vector<int> Nsv,Nsf,Ndv,Ndf;
    igl::circulation(e, true,F,EMAP,EF,EI,Nsv,Nsf);
    igl::circulation(e,false,F,EMAP,EF,EI,Ndv,Ndf);
    old_one_ring = merge(Nsf,Ndf);
  }
  int f1 = EF(e,0);
  int f2 = EF(e,1);
  std::vector<int> new_one_ring = old_one_ring;
  // erase if ==f1 or ==f2
  new_one_ring.erase(
    std::remove(new_one_ring.begin(), new_one_ring.end(), f1), 
    new_one_ring.end());
  new_one_ring.erase(
    std::remove(new_one_ring.begin(), new_one_ring.end(), f2), 
    new_one_ring.end());


  // big box containing new_one_ring
  Eigen::AlignedBox<double,3> big_box;
  // Extend box by placement point
  big_box.extend(p.transpose());
  // Extend box by all other corners (skipping old edge vertices)
  for(const auto f : new_one_ring)
  {
    Eigen::RowVector3d corners[3];
    for(int c = 0;c<3;c++)
    {
      if(F(f,c) == E(e,0) || F(f,c) == E(e,1))
      {
        corners[c] = p;
      }else
      {
        corners[c] = V.row(F(f,c));
        big_box.extend(V.row(F(f,c)).transpose());
      }
    }
    // Degenerate triangles are considered intersections
    if((corners[0]-corners[1]).cross(corners[0]-corners[2]).squaredNorm() < 1e-16)
    {
      return true;
    }
  }
  

  std::vector<const igl::AABB<Eigen::MatrixXd,3>*> candidates;
  tree.append_intersecting_leaves(big_box,candidates);
  

  // Exclude any candidates that are in old_one_ring.
# warning "consider using unordered_set above so that this is O(n+m) rather than O(nm)"
  candidates.erase(
    std::remove_if(candidates.begin(), candidates.end(),
        [&](const igl::AABB<Eigen::MatrixXd,3>* candidate) {
            return std::find(old_one_ring.begin(), old_one_ring.end(), candidate->m_primitive) != old_one_ring.end();
        }),
    candidates.end());
  // print candidates
  constexpr bool stinker = false;
  if(stinker)
  {
    igl::writePLY("before.ply",V,F);
    std::cout<<"Ee = ["<<E(e,0)<<" "<<E(e,1)<<"]+1;"<<std::endl;
    std::cout<<"p = ["<<p<<"];"<<std::endl;
    // print new_one_ring as matlab vector of indices
    std::cout<<"new_one_ring = [";
    for(const auto f : new_one_ring)
    {
      std::cout<<f<<" ";
    }
    std::cout<<"]+1;"<<std::endl;
    // print candidates as matlab vector of indices
    std::cout<<"candidates = [";
    for(const auto * candidate : candidates)
    {
      std::cout<<candidate->m_primitive<<" ";
    }
    std::cout<<"]+1;"<<std::endl;
  }
  
  // For each pair of candidate and new_one_ring, check if they intersect
  bool found_intersection = false;
  for(const int & f : new_one_ring)
  {
    Eigen::AlignedBox<double,3> small_box;
    small_box.extend(p.transpose());
    for(int c = 0;c<3;c++)
    {
      if(F(f,c) != E(e,0) && F(f,c) != E(e,1))
      {
        small_box.extend(V.row(F(f,c)).transpose());
      }
    }
    for(const auto * candidate : candidates)
    {
      const int g = candidate->m_primitive;
      if(!small_box.intersects(candidate->m_box))
      {
        continue;
      }
      // Corner replaced by p
      int c;
      for(c = 0;c<3;c++)
      {
        if(F(f,c) == E(e,0) || F(f,c) == E(e,1))
        {
          break;
        }
      }
      assert(c<3);
      found_intersection = triangle_triangle_intersect(V,F,E,EMAP,EF,EI,f,c,p,g);
      if(found_intersection) { break; }
    }
    if(found_intersection) { break; }
  }
  return found_intersection;
}

///
/// @param[in] orig_pre_collapse  Original pre-collapse callback
/// @param[in] orig_post_collapse Original post-collapse callback
///
#warning "should this take a pad amount?"
void intersection_blocking_collapse_edge_callbacks(
  const igl::decimate_pre_collapse_callback  & orig_pre_collapse,
  const igl::decimate_post_collapse_callback & orig_post_collapse,
        igl::AABB<Eigen::MatrixXd,3> * & tree,
        std::vector<igl::AABB<Eigen::MatrixXd,3> *> & leaves,
        igl::decimate_pre_collapse_callback  & pre_collapse,
        igl::decimate_post_collapse_callback & post_collapse
        )
{
  const double pad = 0;
#warning "is there any reason to capture orig* by _reference_? IIUC capturing by value means we can destroy orig* after this function returns"
  pre_collapse = 
    [orig_pre_collapse,&tree,&leaves](
      const Eigen::MatrixXd & V,
      const Eigen::MatrixXi & F,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & EMAP,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      const igl::min_heap< std::tuple<double,int,int> > & Q,
      const Eigen::VectorXi & EQ,
      const Eigen::MatrixXd & C,
      const int e)->bool
    {
      if(!orig_pre_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e))
      {
        return false;
      }
      
      // Check if there would be (new) intersections
      return 
        !collapse_edge_would_create_intersections(
          e,C.row(e).eval(),V,F,E,EMAP,EF,EI,*tree);
    };
  post_collapse =
    [orig_post_collapse,&tree,&leaves,pad](
      const Eigen::MatrixXd & V,
      const Eigen::MatrixXi & F,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & EMAP,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      const igl::min_heap< std::tuple<double,int,int> > & Q,
      const Eigen::VectorXi & EQ,
      const Eigen::MatrixXd & C,
      const int e,
      const int e1,
      const int e2,
      const int f1,
      const int f2,
      const bool collapsed)
    {
      if(collapsed)
      {
        // detach leaves of deleted faces
        for(int f : {f1,f2})
        {
          auto * sibling = leaves[f]->detach();
          sibling->refit_lineage();
          tree = sibling->root();
          delete leaves[f];
        }
        // If finding `Nf` becomes a bottleneck we could remember it via
        // `pre_collapse` the same way that
        // `qslim_optimal_collapse_edge_callbacks` remembers `v1` and `v2`
        const int m = F.rows();
        const auto survivors = 
          [&m,&e,&EF,&EI,&EMAP](const int f1, const int e1, int & d1)
        {
          int c;
          for(c=0;c<3;c++)
          {
            d1 = EMAP(f1+c*m);
            if((d1 != e) && (d1 != e1)) { break; }
          }
          assert(c<3);
        };
        int d1,d2;
        survivors(f1,e1,d1);
        survivors(f2,e2,d2);
        // Will circulating by continuing in the CCW direction of E(d1,:)
        // encircle the common edge? That is, is E(d1,1) the common vertex?
        const bool ccw = E(d1,1) == E(d2,0) || E(d1,1) == E(d2,1);
#ifndef NDEBUG
        // Common vertex.
        const int s = E(d1,ccw?1:0);
        assert(s == E(d2,0) || s == E(d2,1));
#endif
        std::vector<int> Nf;
        {
          std::vector<int> Nv;
          igl::circulation(d1,ccw,F,EMAP,EF,EI,Nv,Nf);
        }
        for(const int & f : Nf)
        {
          Eigen::AlignedBox<double,3> box;
          box
            .extend(V.row(F(f,0)).transpose())
            .extend(V.row(F(f,1)).transpose())
            .extend(V.row(F(f,2)).transpose());
          // Always grab root (returns self if no update)
          tree = leaves[f]->update(box,pad)->root();
          assert(tree == tree->root());
        }
          assert(tree == tree->root());
      }
//#ifndef NDEBUG
#if false
#warning "🐌🐌🐌🐌🐌🐌🐌🐌🐌🐌🐌 Slow intersection checking..."
      constexpr bool stinker = true;
      if(stinker && igl::copyleft::cgal::is_self_intersecting(V,F))
      {
        igl::writePLY("after.ply",V,F);
        printf("💩💩💩💩💩 Just shit the bed on e=%d 🛌🛌🛌🛌 \n",e);
        printf("🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥\n");
        printf("🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥\n");
        printf("🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥\n");
        printf("🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨\n");
        printf("🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨\n");
        printf("🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨🧨\n");
        printf("💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥\n");
        printf("💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥\n");
        printf("💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥💥\n");
        printf("☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ \n");
        printf("☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ \n");
        printf("☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ ☠️ \n");
        exit(1);
      }
#endif
      // Finally. Run callback.
      return orig_post_collapse(
        V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2,collapsed);
    };
}
void intersection_blocking_collapse_edge_callbacks(
  igl::AABB<Eigen::MatrixXd,3> * & tree,
  std::vector<igl::AABB<Eigen::MatrixXd,3> *> & leaves,
  igl::decimate_pre_collapse_callback  & pre_collapse,
  igl::decimate_post_collapse_callback & post_collapse)
{
  igl::decimate_pre_collapse_callback  always_try;
  igl::decimate_post_collapse_callback never_care;
  igl::decimate_trivial_callbacks(always_try,never_care);
  intersection_blocking_collapse_edge_callbacks(
    always_try,
    never_care,
    tree,
    leaves,
    pre_collapse,
    post_collapse);
}

int main(int argc, char *argv[])
{
  IGL_TICTOC_LAMBDA;
#if true
  Eigen::MatrixXd V,V0;
  Eigen::MatrixXi F,F0;
  const bool use_test = argc<=1;
  igl::read_triangle_mesh(use_test?"/Users/alecjacobson/Downloads/simple-intersection-collapse-2.ply":argv[1],V,F);
  V0 = V;F0 = F;
  ///////////////////////////////////////////////////////////
  /// Before collapsing starts
  ///////////////////////////////////////////////////////////
  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  printf("tree->init: %g\n",tictoc());
  validate(tree);
#warning "should pad leaves?"
  // Gather list of pointers to leaves
  std::vector<igl::AABB<Eigen::MatrixXd,3>*> leaves = tree->gather_leaves(F.rows());

  // Dummy
  if(use_test)
  {
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E,EF,EI;
    igl::edge_flaps(F,E,EMAP,EF,EI);
    igl::decimate_pre_collapse_callback  pre_collapse;
    igl::decimate_post_collapse_callback post_collapse;
    intersection_blocking_collapse_edge_callbacks(
      tree,
      leaves,
      pre_collapse,
      post_collapse);
    
    igl::min_heap< std::tuple<double,int,int> > Q;
    Eigen::VectorXi EQ;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(),3);

    // Try to do the collapses.
    std::vector<int> edges_to_collapse;
    for(const auto edge : {Eigen::RowVector2i(3,6),Eigen::RowVector2i(6,18)})
    {
      int e = -1;
      const bool found = 
        ( (E.array().col(0)==edge[0] && E.array().col(1)==edge[1])||
          (E.array().col(0)==edge[1] && E.array().col(1)==edge[0])).maxCoeff(&e);
      assert(found);
      edges_to_collapse.push_back(e);
    }
    for(auto & e : edges_to_collapse)
    {
      printf("E(e=%d,:) = %d,%d\n",e,E(e,0),E(e,1));
      bool collapsed = true;
      int e1,e2,f1,f2;
      // pre_collapse only has access to C not p
      C.row(e) = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
      if(pre_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e))
      {
        collapsed = igl::collapse_edge(
          e,C.row(e).eval(),
          /*Nsv,Nsf,Ndv,Ndf,*/
          V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
      }else
      {
        // Aborted by pre collapse callback
        collapsed = false;
      }
      printf("!!!!!!!!!!!!!!!!!collapsed? %s\n",collapsed?"✅":"❌");
      post_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2,collapsed);
    }
    printf("tree: %p\n",tree);
  }else
  {
    const int target_m = F.rows() * 0.1 + 1;
    for(auto pass : {0,1})
    {
      igl::decimate_pre_collapse_callback  pre_collapse;
      igl::decimate_post_collapse_callback post_collapse;
      igl::decimate_trivial_callbacks(pre_collapse,post_collapse);
      if(pass == 1)
      {
        intersection_blocking_collapse_edge_callbacks(
          pre_collapse, post_collapse, // These will get copied as needed
          tree, leaves,
          pre_collapse, post_collapse);
      }
      Eigen::MatrixXd VO = V;
      Eigen::MatrixXi FO = F;
#     warning "Should connect_boundary_to_infinity"
      Eigen::VectorXi EMAP;
      Eigen::MatrixXi E,EF,EI;
      igl::edge_flaps(FO,E,EMAP,EF,EI);
      int m = F.rows();
      const int orig_m = m;
      Eigen::MatrixXd U;
      Eigen::MatrixXi G;
      Eigen::VectorXi J,I;
      tictoc();
      const bool ret = igl::decimate(
        VO, FO,
        igl::shortest_edge_and_midpoint,
        igl::max_faces_stopping_condition(m,orig_m,target_m),
        pre_collapse,
        post_collapse,
        E, EMAP, EF, EI,
        U, G, J, I);
      printf("pass %d in %g sec\n",pass,tictoc());
      igl::writePLY(STR("out-"<<pass<<".ply"),U,G);
      if(pass == 1)
      {
        V = U;
        F = G;
      }
    }
  }

  vis(V,F,*tree);

  validate(tree);
  tree = tree->root();
  delete tree;
  leaves.clear();
  exit(1);

  //Eigen::MatrixXd K = Eigen::MatrixXd::Constant(F.rows(),3,0.9);
  //for(const auto f : new_one_ring)
  //{
  //  K.row(f) = Eigen::RowVector3d(0.3,0.9,0.3);
  //}
  //Eigen::MatrixXd K0 = Eigen::MatrixXd::Constant(F.rows(),3,0.9);
  //for(const auto f : old_one_ring)
  //{
  //  K0.row(f) = Eigen::RowVector3d(0.9,0.3,0.3);
  //}


  igl::opengl::glfw::Viewer vr;
  //vr.data().set_mesh(V0,F0);
  //vr.data().set_face_based(true);
  //vr.data().set_colors(K0);
  //vr.data().show_faces = false;
  //vr.data().line_color = Eigen::RowVector4f(1,1,1,1);
  //vr.data().double_sided = true;
  //vr.append_mesh();
  vr.data().set_mesh(V,F);
  //vr.data().set_colors(K);
  vr.data().set_face_based(true);
  vr.data().show_lines = false;
  vr.data().double_sided = true;

  //Eigen::MatrixXd TV;
  //Eigen::MatrixXi TQ,TE;
  //box_faces(big_box,0,TV,TQ);
  //quad_edges(TQ,TE);
  //vr.append_mesh();
  //vr.data().set_edges(TV,TE,Eigen::RowVector3d(1,1,1));
  //vr.data().line_width = 2;

  
  vr.launch();

#else
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc>1?argv[1]:TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  // make into soup
  V = V(Eigen::Map<Eigen::VectorXi>(F.data(),F.size()), Eigen::all).eval();
  F = Eigen::Map<Eigen::MatrixXi>(igl::colon<int>(0,V.rows()-1).data(),V.rows()/3,3).eval();

  //F = F.topRows(4).eval();
  //
  Eigen::MatrixXd BC;
  igl::barycenter(V,F,BC);
  Eigen::VectorXd sqrD;
  Eigen::VectorXi I;
  Eigen::VectorXi J = igl::colon<int>(0,F.rows()-1);
  Eigen::MatrixXd C;

  IGL_TICTOC_LAMBDA;
  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  printf("tree.init(): %g secs\n",tictoc());
  tictoc();
  tree->squared_distance(V,F,BC,sqrD,I,C);
  assert(I.isApprox(J,0));
  printf("tree.squared_distance(): %g secs\n",tictoc());
  printf("  surface_area: %g\n",tree->internal_surface_area());
  printf("  is_root(): %d\n",tree->is_root());
  printf("  size: %d\n",size(tree));
  printf("  height: %d/%d\n",height(tree),(int)F.rows());
  //print(&tree);
  validate(tree);

  // Gather list of pointers to leaves
  std::vector<igl::AABB<Eigen::MatrixXd,3>*> leaves = tree->gather_leaves(F.rows());
  for(auto * leaf : leaves) { assert(leaf); }
  //print(tree);
  printf("--------------------------------\n");

  // detach and insert each leaf
  const double h = igl::avg_edge_length(V,F);
  const double pad = h;
  {
    // Gather list of pointers to leaves
    std::vector<igl::AABB<Eigen::MatrixXd,3>*> leaves = tree->gather_leaves(F.rows());
    for(auto * leaf : leaves) { assert(leaf); }
    tree = tree->pad(leaves,pad,2);
  }

  tictoc();
  tree->squared_distance(V,F,BC,sqrD,I,C);
  assert(I.isApprox(J,0));
  printf("tree.squared_distance(): %g secs\n",tictoc());
  printf("  surface_area: %g\n",tree->internal_surface_area());
  printf("  is_root(): %d\n",tree->is_root());
  printf("  size: %d\n",size(tree));
  printf("  height: %d/%d\n",height(tree),(int)F.rows());
  //print(tree);


  Eigen::VectorXi RV;
  Eigen::VectorXi RF;
  // Perturb a small subset of the triangles
  {
    {
      igl::randperm(F.rows(),RF);
      RF = RF.topRows(std::min(12,(int)F.rows())).eval();
      RV.resize(RF.size()*3);
      RV << RF, RF.array()+F.rows(), RF.array()+2*F.rows();
    }
    Eigen::MatrixXd TF = 0.1*h*Eigen::MatrixXd::Random(RF.size(),3);
    Eigen::MatrixXd TV(RV.rows(),3);
    TV<<TF,TF,TF;
    V(RV,Eigen::all) += TV;
    igl::barycenter(V,F,BC);
  }
  const int qi = RF(0);


  {
    tictoc();
    for(int i = 0;i<RF.size();i++)
    {
      tree = leaves[RF(i)]->update_primitive(V,F,pad)->root();
    }
    printf("        tree.refit                :          %g secs\n",tictoc());
    tictoc();
    for(int r = 0;r<MAX_RUNS;r++)
    {
      tree->squared_distance(V,F,Eigen::MatrixXd(BC.row(qi)),sqrD,I,C);
    }
    assert(I(0) == qi);
    printf("%d,%g ← tree.squared_distance(0,…):          %g secs\n",I(0),sqrD(0),tictoc()/MAX_RUNS);
  }
  {
    tictoc();
    for(int r = 0;r<MAX_RUNS;r++)
    {
      igl::point_mesh_squared_distance(Eigen::MatrixXd(BC.row(qi)),V,F,sqrD,I,C);
    }
    assert(I(0) == qi);
    printf("%d,%g ← point_mesh_squared_distance(0,…):    %g secs\n",I(0),sqrD(0),tictoc()/MAX_RUNS);
  }

  {
    tictoc();
    double min_dist;
    int min_i;
    for(int r = 0;r<MAX_RUNS;r++)
    {
      min_dist = std::numeric_limits<double>::infinity();
      min_i = -1;
      Eigen::RowVector3d q = Eigen::RowVector3d(BC.row(qi));
      for(int i = 0;i<F.rows();i++)
      {
        Eigen::RowVector3d c;
        double d;
        igl::point_simplex_squared_distance<3>(q,V,F,i,d,c);
        if(d < min_dist)
        {
          min_dist = d;
          min_i = i;
        }
      }
    }
    printf("%d,%g ← point_simplex_squared_distance(0,…): %g secs\n",min_i,min_dist,tictoc()/MAX_RUNS);
    assert(min_i == qi);
  }

  
  vis(V,F,*tree);
  delete tree;


  //igl::AABB<Eigen::MatrixXd, 3> * dynamic = nullptr;
  //std::vector<igl::AABB<Eigen::MatrixXd, 3> > leaves;
  //for(auto rotation_amount : {1, 2, 3})
  //{
  //  // if vector is storing objects, must clear first.
  //  leaves.clear();
  //  // tree is now invalid, but deleting is safe.
  //  delete dynamic;
  //  {
  //    printf("\n--------------------------------\n\n");
  //    printf("rotation_amount: %d\n",rotation_amount);
  //    tictoc();
  //    // The root starts as the first one which will be self-inserted
  //    leaves.resize(F.rows());
  //    dynamic = leaves.data();
  //    for(int i = 0;i<F.rows();i++)
  //    {
  //      auto * leaf = &leaves[i];
  //      // Use the idiotic .init()
  //      leaf->init(V,F,Eigen::MatrixXi(),(Eigen::VectorXi(1)<<i).finished());
  //      dynamic = dynamic->insert(leaf)->root();

  //      if(rotation_amount==1)
  //      {
  //        const bool ret = leaf->rotate();
  //      }

  //      if(rotation_amount>=2)
  //      {
  //        leaf->rotate_lineage();
  //      }
  //    }
  //    if(rotation_amount==3)
  //    {
  //      std::deque<igl::AABB<Eigen::MatrixXd, 3> *> bfs;
  //      bfs.push_back(dynamic);
  //      std::vector<igl::AABB<Eigen::MatrixXd, 3> *> all_nodes;
  //      while(!bfs.empty())
  //      {
  //        auto * node = bfs.back();
  //        bfs.pop_back();
  //        if(node->m_left)
  //        {
  //          bfs.push_back(node->m_left);
  //        }
  //        if(node->m_right)
  //        {
  //          bfs.push_back(node->m_right);
  //        }
  //        all_nodes.push_back(node);
  //      }
  //      while(!all_nodes.empty())
  //      {
  //        auto * node = all_nodes.back();
  //        all_nodes.pop_back();
  //        assert(node);
  //        const bool ret = node->rotate();
  //      }
  //    }
  //    printf("dynamic %g\n",tictoc());
  //    tictoc();
  //    dynamic->squared_distance(V,F,BC,sqrD,I,C);
  //    assert(I.isApprox(J,0));
  //    printf("tree.squared_distance(): %g secs\n",tictoc());
  //    printf("  surface_area: %g\n",dynamic->internal_surface_area());
  //    printf("  is_root(): %d\n",dynamic->is_root());
  //    printf("  size: %d\n",size(dynamic));
  //    printf("  height: %d/%d\n",height(dynamic),(int)F.rows());
  //    //print(dynamic);
  //    validate(dynamic, leaves);
  //  }
  //  printf("********\n");
  //}

  //vis(V,F,*dynamic);

#endif
}
