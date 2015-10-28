#ifndef IGL_CGAL_EXTRACT_CELLS
#define IGL_CGAL_EXTRACT_CELLS

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
    namespace cgal {
        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedP,
            typename DeriveduE,
            typename uE2EType,
            typename DerivedEMAP,
            typename DerivedC >
        IGL_INLINE size_t extract_cells_single_component(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedP>& P,
                const Eigen::PlainObjectBase<DeriveduE>& uE,
                const std::vector<std::vector<uE2EType> >& uE2E,
                const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
                Eigen::PlainObjectBase<DerivedC>& cells);

        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedP,
            typename DerivedE,
            typename DeriveduE,
            typename uE2EType,
            typename DerivedEMAP,
            typename DerivedC >
        IGL_INLINE size_t extract_cells(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedP>& P,
                const Eigen::PlainObjectBase<DerivedE>& E,
                const Eigen::PlainObjectBase<DeriveduE>& uE,
                const std::vector<std::vector<uE2EType> >& uE2E,
                const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
                Eigen::PlainObjectBase<DerivedC>& cells);

        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedC >
        IGL_INLINE size_t extract_cells(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                Eigen::PlainObjectBase<DerivedC>& cells);
    }
}

#ifndef IGL_STATIC_LIBRARY
#  include "extract_cells.cpp"
#endif
#endif
