#include <igl/polyvector_field_cut_mesh_with_singularities.h>
#include <igl/dijkstra.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/is_border_vertex.h>
#include <igl/cut_mesh_from_singularities.h>
#include <set>


template <typename DerivedV, typename DerivedF, typename VFType, typename VVType, typename DerivedTT, typename DerivedC, typename DerivedS>
IGL_INLINE void igl::polyvector_field_cut_mesh_with_singularities(
                                                                  const Eigen::PlainObjectBase<DerivedV> &V,
                                                                  const Eigen::PlainObjectBase<DerivedF> &F,
                                                                  const std::vector<std::vector<VFType> >& VF,
                                                                  const std::vector<std::vector<VVType> >& VV,
                                                                  const Eigen::PlainObjectBase<DerivedTT>& TT,
                                                                  const Eigen::PlainObjectBase<DerivedTT>& TTi,
                                                                  const Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                  Eigen::PlainObjectBase<DerivedC> &cuts)
{

  //first, get a spanning tree for the mesh (no missmatch needed)
  igl::cut_mesh_from_singularities(V, F, Eigen::MatrixXd::Zero(F.rows(), 3).eval(), cuts);

  std::set<int> vertices_in_cut;
  for (int i =0; i< cuts.rows(); ++i)
    for (int j =0;j< cuts.cols(); ++j)
      if (cuts(i,j))
        vertices_in_cut.insert(F(i,j));

  //then, add all singularities one by one by using Dijkstra's algorithm
  for (int i = 0; i<singularities.rows(); ++i)
  {
    std::vector<int> path;
    Eigen::VectorXd min_distance;
    Eigen::VectorXi previous;
    int vertex_found = igl::dijkstra_compute_paths(singularities[i], vertices_in_cut, VV, min_distance, previous);
    if(vertex_found ==-1)
      // this means that there are no cuts
      path.push_back(singularities[i]);
    else
      igl::dijkstra_get_shortest_path_to(vertex_found, previous, path);

    vertices_in_cut.insert(path.begin(), path.end());

    //insert to cut
    for (int ii = 0; ii<path.size()-1; ++ii)
    {
      const int &v0 = path[ii];
      const int &v1 = path[ii+1];

      std::vector<int> vf0 = VF[v0];
      std::sort(vf0.begin(), vf0.end());
      std::vector<int> vf1 = VF[v1];
      std::sort(vf1.begin(), vf1.end());
      std::vector<int> common_face_v(std::max(vf0.size(),vf1.size()));
      std::vector<int>::iterator it;
      it=std::set_intersection (vf0.begin(), vf0.end(), vf1.begin(), vf1.end(), common_face_v.begin());
      common_face_v.resize(it-common_face_v.begin());
      assert(common_face_v.size() == 2);

      const int &fi = common_face_v[0];
      int j=-1;
      for (unsigned z=0; z<3; ++z)
        if (((F(fi,z) == v0) && (F(fi,(z+1)%3) == v1)) ||((F(fi,z) == v1) && (F(fi,(z+1)%3) == v0)))
          j=z;
      assert(j!=-1);
      cuts(fi,j) = 1;
      cuts(TT(fi,j), TTi(fi,j)) = 1;

    }
  }

}


//Wrapper of the above with only vertices and faces as mesh input
template <typename DerivedV, typename DerivedF, typename DerivedC, typename DerivedS>
IGL_INLINE void igl::polyvector_field_cut_mesh_with_singularities(
                                                                  const Eigen::PlainObjectBase<DerivedV> &V,
                                                                  const Eigen::PlainObjectBase<DerivedF> &F,
                                                                  const Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                  Eigen::PlainObjectBase<DerivedC> &cuts)
{

  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  std::vector<std::vector<int> > VV;
  igl::adjacency_list(F, VV);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(V,F,TT,TTi);

  igl::polyvector_field_cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities, cuts);


}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::polyvector_field_cut_mesh_with_singularities<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, int, int, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::polyvector_field_cut_mesh_with_singularities<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
