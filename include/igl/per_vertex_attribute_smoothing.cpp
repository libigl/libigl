#include "per_vertex_attribute_smoothing.h"
#include <vector>

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::per_vertex_attribute_smoothing(
    const Eigen::PlainObjectBase<DerivedV>& Ain,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV> & Aout)
{
    std::vector<double> denominator(Ain.rows(), 0);
    Aout = Eigen::PlainObjectBase<DerivedV>::Zero(Ain.rows(), Ain.cols());
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int j1 = (j + 1) % 3;
            int j2 = (j + 2) % 3;
            Aout.row(F(i, j)) += Ain.row(F(i, j1)) + Ain.row(F(i, j2));
            denominator[F(i, j)] += 2;
        }
    }
    for (int i = 0; i < Ain.rows(); ++i)
        Aout.row(i) /= denominator[i];
}

