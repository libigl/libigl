#ifndef ORDER_FACETS_AROUND_EDGES
#define ORDER_FACETS_AROUND_EDGES

#include "../igl_inline.h"
#include <Eigen/Core>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace igl {
    // For each undirected edge, sort its adjacent faces.
    //
    // Inputs:
    //   V    #V by 3 list of vertices.
    //   F    #F by 3 list of faces
    //   N    #F by 3 list of face normals.
    //   E    #F*3 by 2 list vertex indices, represents directed edges.
    //  uE    #uE by 2 list of vertex_indices, represents undirected edges.
    //  EMAP  #F*3 list of indices that maps E to uE. (a many-to-one map)
    //  uE2E  #uE list of lists that maps uE to E. (a one-to-many map)
    //
    // Outputs:
    //   uE2oE #uE list of lists that maps uE to an ordered list of E. (a
    //         one-to-many map)
    //   uE2C  #uE list of lists of bools indicates whether each face in
    //         uE2oE[i] is consistently orientated as the ordering.
    //
    template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedN,
        typename DerivedE,
        typename DeriveduE,
        typename DerivedEMAP,
        typename uE2EType,
        typename uE2oEType,
        typename uE2CType >
    IGL_INLINE
    typename std::enable_if<!std::is_same<typename DerivedV::Scalar,
    typename CGAL::Exact_predicates_exact_constructions_kernel::FT>::value, void>::type
    order_facets_around_edges(
            const Eigen::PlainObjectBase<DerivedV>& V,
            const Eigen::PlainObjectBase<DerivedF>& F,
            const Eigen::PlainObjectBase<DerivedN>& N,
            const Eigen::PlainObjectBase<DerivedE>& E,
            const Eigen::PlainObjectBase<DeriveduE>& uE,
            const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
            const std::vector<std::vector<uE2EType> >& uE2E,
            std::vector<std::vector<uE2oEType> >& uE2oE,
            std::vector<std::vector<uE2CType > >& uE2C );

    template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedN,
        typename DerivedE,
        typename DeriveduE,
        typename DerivedEMAP,
        typename uE2EType,
        typename uE2oEType,
        typename uE2CType >
    IGL_INLINE 
    typename std::enable_if<std::is_same<typename DerivedV::Scalar,
    typename CGAL::Exact_predicates_exact_constructions_kernel::FT>::value, void>::type
    order_facets_around_edges(
            const Eigen::PlainObjectBase<DerivedV>& V,
            const Eigen::PlainObjectBase<DerivedF>& F,
            const Eigen::PlainObjectBase<DerivedN>& N,
            const Eigen::PlainObjectBase<DerivedE>& E,
            const Eigen::PlainObjectBase<DeriveduE>& uE,
            const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
            const std::vector<std::vector<uE2EType> >& uE2E,
            std::vector<std::vector<uE2oEType> >& uE2oE,
            std::vector<std::vector<uE2CType > >& uE2C );
}

#ifndef IGL_STATIC_LIBRARY
#include "order_facets_around_edges.cpp"
#endif

#endif
