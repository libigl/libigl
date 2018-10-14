#include "generate_sparse_marching_cubes_grid.h"

#include <unordered_map>
#include <array>
#include <vector>


struct RowVector3iHash {
    std::size_t operator()(const Eigen::RowVector3i& key) const {
        std::size_t seed = 0;
        for (int i = 0; i < 3; i++) {
            seed ^= key[i] + 0x9e3779b9 + (seed<<6) + (seed>>2); // Copied from boost::hash_combine
        }
        return std::hash<std::size_t>()(seed);
    }
};

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template<class DerivedM, class DerivedR>
int matrix_push_back(Eigen::PlainObjectBase<DerivedM>& mat, const Eigen::PlainObjectBase<DerivedR>& row, int n) {
    if (mat.rows() == n) {
        mat.conservativeResize(2*n+1, mat.cols());
    }

    mat.row(n) = row;
    return n + 1;
}


template<class DerivedM>
int vector_push_back(Eigen::PlainObjectBase<DerivedM>& mat, typename Eigen::PlainObjectBase<DerivedM>::Scalar row, int n) {
    if (mat.rows() == n) {
        mat.conservativeResize(2*n+1);
    }

    mat[n] = row;
    return n + 1;
}



template <typename DerivedS, typename DerivedV, typename DerivedI>
void igl::generate_sparse_marching_cubes_grid(const Eigen::RowVector3d& p0, std::function<double(const Eigen::RowVector3d&)> scalarFunc, double eps,
                                              Eigen::PlainObjectBase<DerivedS>& CS, Eigen::PlainObjectBase<DerivedV>& CV, Eigen::PlainObjectBase<DerivedI>& CI) {
    using namespace std;
    using namespace Eigen;

    double half_eps = 0.5 * eps;
    int num_verts = 0, num_cubes = 0;
    CV.resize(0, 3);
    CS.resize(0, 1);
    CI.resize(0, 8);

    unordered_map<RowVector3i, int, RowVector3iHash> visited;
    vector<RowVector3i> queue { RowVector3i(0, 0, 0) };

    while (queue.size() > 0) {
        RowVector3i pi = queue.back();
        queue.pop_back();

        RowVector3d ctr = p0 + eps*pi.cast<double>(); // R^3 center of this cube

        // X, Y, Z basis vectors, and array of neighbor offsets
        const RowVector3i bx(1, 0, 0), by(0, 1, 0), bz(0, 0, -1); // I flipped the z-axis to because I'm a derp and got the face winding backwards :)
        const array<RowVector3i, 6> neighbors = {
            bx, -bx, by, -by, bz, -bz
        };

        // Compute the position of the cube corners and the scalar values at those corners
        array<RowVector3d, 8> cubeCorners = {
          ctr+half_eps*(bx+by+bz).cast<double>(), ctr+half_eps*(bx+by-bz).cast<double>(), ctr+half_eps*(-bx+by-bz).cast<double>(), ctr+half_eps*(-bx+by+bz).cast<double>(),
          ctr+half_eps*(bx-by+bz).cast<double>(), ctr+half_eps*(bx-by-bz).cast<double>(), ctr+half_eps*(-bx-by-bz).cast<double>(), ctr+half_eps*(-bx-by+bz).cast<double>()
        };
        array<double, 8> cubeScalars;
        for (int i = 0; i < 8; i++) { cubeScalars[i] = scalarFunc(cubeCorners[i]); }

        // If this cube doesn't intersect the surface, disregard it
        bool validCube = false;
        double sign = sgn(cubeScalars[0]);
        for (int i = 1; i < 8; i++) {
            if (sign != sgn(cubeScalars[i])) {
                validCube = true;
                break;
            }
        }
        if (!validCube) {
            continue;
        }


        // Add the cube vertices and indices to the output arrays if they are not there already
        RowVectorXi cube(1, 8);
        uint8_t vertexAlreadyAdded = 0;
        constexpr array<uint8_t, 6> zv = {
            (1 << 0) | (1 << 1) | (1 << 4) | (1 << 5),
            (1 << 2) | (1 << 3) | (1 << 6) | (1 << 7),
            (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3),
            (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7),
            (1 << 0) | (1 << 3) | (1 << 4) | (1 << 7),
            (1 << 1) | (1 << 2) | (1 << 5) | (1 << 6),
        };
        constexpr array<array<int, 4>, 6> zvv {{
            {{0, 1, 4, 5}}, {{3, 2, 7, 6}}, {{0, 1, 2, 3}},
            {{4, 5, 6, 7}}, {{0, 3, 4, 7}}, {{1, 2, 5, 6}}
        }};

        for (int n = 0; n < 6; n++) { // For each neighbor, check if its added before
            RowVector3i nkey = pi + neighbors[n];
            auto nbr = visited.find(nkey);
            if (nbr != visited.end()) { // We've already visited this neighbor
                vertexAlreadyAdded |= zv[n];
                for (int i = 0; i < 4; i++) { cube[zvv[n][i]] = CI.row(nbr->second)[zvv[n % 2 == 0 ? n + 1 : n - 1][i]]; }
            } else {
                queue.push_back(nkey);
            }
        }

        for (int i = 0; i < 8; i++) { // Add new vertices to the arrays
            if (0 == ((1 << i) & vertexAlreadyAdded)) {
                cube[i] = num_verts;
                matrix_push_back(CV, cubeCorners[i], num_verts);
                num_verts = vector_push_back(CS, cubeScalars[i], num_verts);
            }
        }

        visited[pi] = num_cubes;
        num_cubes = matrix_push_back(CI, cube, num_cubes);
    }

    CV.conservativeResize(num_verts, 3);
    CS.conservativeResize(num_verts);
    CI.conservativeResize(num_cubes, 8);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template
void igl::generate_sparse_marching_cubes_grid<
    Eigen::Matrix<double, -1, 1, 0, -1, 1>,
    Eigen::Matrix<double, -1, -1, 0, -1, -1>,
    Eigen::Matrix<int, -1, -1, 0, -1, -1>>
        (const Eigen::RowVector3d& p0,  std::function<double(const Eigen::RowVector3d&)> scalarFunc, double eps,
         Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>& CS,
         Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>& CV,
         Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>& CI);
#endif
