#pragma once
#include <string>

#include <Eigen/Core>
#include <gtest/gtest.h>

#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>

namespace test_common {
    template<typename DerivedV, typename DerivedF>
    void load_mesh(
            const std::string& filename, 
            Eigen::PlainObjectBase<DerivedV>& V,
            Eigen::PlainObjectBase<DerivedF>& F) {
        auto find_file = [&](const std::string& val) {
            return std::string(TEST_DIR) + "/data/" + val;
        };
        igl::read_triangle_mesh(find_file(filename), V, F);
    }

    template<typename Derived>
    void load_matrix(const std::string& filename,
            Eigen::PlainObjectBase<Derived>& M) {
        igl::readDMAT(std::string(TEST_DIR) + "/data/" + filename, M);
    }
}
