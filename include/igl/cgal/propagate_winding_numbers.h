#ifndef IGL_CGAL_PROPAGATE_WINDING_NUMBERS_H
#define IGL_CGAL_PROPAGATE_WINDING_NUMBERS_H
#include "../igl_inline.h"
#include <Eigen/Core>
#include <vector>

// The following methods compute the winding number on each side of each facet
// or patch of a 3D mesh.  The input mesh is valid if it splits the ambient
// space, R^3, into subspaces with constant integer winding numbers.  That is
// the input mesh should be closed and for each directed edge the number of
// clockwise facing facets should equal the number of counterclockwise facing
// facets.

namespace igl {
    namespace cgal {
        // Compute winding number on each side of each patch.  All patches must
        // be connected.  The per-patch label partitions the patches into
        // different components, and winding number is computed for each
        // component.  For example, ``patch_W(i, j*2)`` for ``i`` in [0, #P)
        // and ``j`` in [0, k) is the winding number on the positive side of
        // patch ``i`` with respect to component ``j``.  Similarly,
        // ``patch_W(i, j*2+1)`` is the winding number on the negative side of
        // patch ``i`` with respect to component ``j``.
        //
        // Patch and intersection curve can be generated using:
        // ``igl::extract_manifold_patches(...)`` and
        // ``igl::extract_non_manifold_edge_curves(...)``
        //
        //
        // Inputs:
        //   V  #V by 3 list of vertices.
        //   F  #F by 3 list of trianlges.
        //   uE #uE by 2 list of undirected edges.
        //   uE2E  #uE list of lists mapping each undirected edge to directed
        //         edges.
        //   labels  #P list of labels, ranging from 0 to k-1.
        //           These labels partition patches into k components, and
        //           winding number is computed for each component.
        //   P  #F list of patch indices.
        //   intersection_curves  A list of non-manifold curves that separates
        //                        the mesh into patches.
        //
        // Outputs:
        //   patch_W  #P by k*2 list of winding numbers on each side of a
        //            patch for each component.
        //
        // Returns:
        //   True iff integer winding numbers can be consistently assigned.
        template<
            typename DerivedV,
            typename DerivedF,
            typename DeriveduE,
            typename uE2EType,
            typename DerivedC,
            typename DerivedP,
            typename DerivedW >
        IGL_INLINE bool propagate_winding_numbers_single_component_patch_wise(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DeriveduE>& uE,
                const std::vector<std::vector<uE2EType> >& uE2E,
                const Eigen::PlainObjectBase<DerivedC>& labels,
                const Eigen::PlainObjectBase<DerivedP>& P,
                const std::vector<std::vector<size_t> >& intersection_curves,
                Eigen::PlainObjectBase<DerivedW>& patch_W);

        // Compute winding number on each side of the face.  The input mesh must
        // be a single edge-connected component.
        //
        // Inputs:
        //   V  #V by 3 list of vertex positions.
        //   F  #F by 3 list of triangle indices into V.
        //   labels  #F list of facet labels ranging from 0 to k-1.
        //
        // Output:
        //   W  #F by 2*k list of winding numbers.
        //      ``W(i, 2*j)`` is the winding number on the positive side of
        //      facet ``i`` with respect to facets labelled ``j``.
        //      ``W(i, 2*j+1)`` is the winding number on the negative side of
        //      facet ``i`` with respect to facets labelled ``j``.
        //
        // Returns:
        //   True iff integer winding number can be consistently assigned.
        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedC,
            typename DerivedW>
        IGL_INLINE bool propagate_winding_numbers_single_component(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedC>& labels,
                Eigen::PlainObjectBase<DerivedW>& W);

        // Compute the winding number on each side of the face.
        // This method assumes the input mesh (V, F) forms a single connected
        // component.
        //
        // Inputs:
        //   V  #V by 3 list of vertex positions.
        //   F  #F by 3 list of triangle indices into V.
        //
        // Output:
        //   W  #F by 2 list of winding numbers.  W(i,0) is the winding number
        //   on the positive side of facet i, and W(i, 1) is the winding
        //   number on the negative side of facet I[i].
        //
        // Returns:
        //   True iff integer winding number can be consistently assigned.
        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedW>
        IGL_INLINE bool propagate_winding_numbers_single_component(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                Eigen::PlainObjectBase<DerivedW>& W);

        // Compute winding number on each side of the face.  The input mesh
        // could contain multiple connected components.  The input mesh must
        // represent the boundary of a valid 3D volume, which means it is
        // closed, consistently oriented and induces integer winding numbers.
        //
        // Inputs:
        //   V  #V by 3 list of vertex positions.
        //   F  #F by 3 list of triangle indices into V.
        //   labels  #F list of facet labels ranging from 0 to k-1.
        //
        // Output:
        //   W  #F by k*2 list of winding numbers.  ``W(i,j*2)`` is the winding
        //      number on the positive side of facet ``i`` with respect to the
        //      facets labeled ``j``.  Similarly, ``W(i,j*2+1)`` is the winding
        //      number on the negative side of facet ``i`` with respect to the
        //      facets labeled ``j``.
        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedL,
            typename DerivedW>
        IGL_INLINE void propagate_winding_numbers(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedL>& labels,
                Eigen::PlainObjectBase<DerivedW>& W);
    }
}

#ifndef IGL_STATIC_LIBRARY
#  include "propagate_winding_numbers.cpp"
#endif
#endif
