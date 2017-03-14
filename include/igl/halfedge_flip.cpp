
#include "halfedge_flip.h"
#include "parallel_for.h"
#include "sortrows.h"
#include <iostream>

#define min(x, y) ((x < y) ? x : y)
#define max(x, y) ((x > y) ? x : y) 

template <typename DerivedF, typename DerivedHE>
IGL_INLINE void igl::halfedge_flip(
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedHE>& HE)
{
    using namespace std;
    const int m = F.rows();
    switch(F.cols())
    {
        case 3:
        {
            /* each row is (v_start, v_target, f_id, e_rank) */
            Eigen::MatrixXi H(3 * m, 4);
            parallel_for(
                m,
                [&F, &H](const int i)
                {
                    int v0 = F(i, 0);
                    int v1 = F(i, 1);
                    int v2 = F(i, 2);

                    H.row(3 * i    ) << min(v1, v2), max(v1, v2), i, 0;
                    H.row(3 * i + 1) << min(v2, v0), max(v2, v0), i, 1;
                    H.row(3 * i + 2) << min(v0, v1), max(v0, v1), i, 2;
                },
                1000);

            /* sort H and then get halfedge pairs */
            Eigen::MatrixXi H_sort;
            Eigen::VectorXi I;
            sortrows(H, true, H_sort, I);
            
            /* extract opposite adjacency */
            HE = - Eigen::PlainObjectBase<DerivedHE>::Ones(m, 3);
            int iter = 1;
            while (iter < 3 * m)
            {
                /* if the halfedge matches, write to adjT and adjR */
                if (H_sort(iter, 0) == H_sort(iter - 1, 0) &&
                    H_sort(iter, 1) == H_sort(iter - 1, 1))
                {
                    /* get useful information (f_id, v_rank, v_id) */
                    auto a = H_sort.row(iter).tail(2);
                    auto b = H_sort.row(iter - 1).tail(2);

                    /* record halfedge flip relation */
                    HE(a(0), a(1)) = 3 * b(0) + b(1);
                    HE(b(0), b(1)) = 3 * a(0) + a(1);

                    iter += 2;
                }else{
                    iter++;
                }
            }
            break;
        }
        default:
        {
            cerr << "opposite_adjacency.h: Error: Simplex size ("<< F.cols() <<
                ") not supported" << endl;
            assert(false);
        }
    }
};
