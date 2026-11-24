#pragma once
#include <Eigen/Core>
#include <vector>

using Point = Eigen::RowVector3d;
struct Triangle { Point a, b, c; };
struct Ray { Point origin, direction; };

#define DECLARE_ADAPTER(name) \
namespace name { \
int classify(const Ray &, const Triangle &); \
int compare(const Ray &, const Triangle &, const Triangle &); \
bool first(const Ray &, const Eigen::MatrixXd &, const Eigen::MatrixXi &, int &); \
bool all(const Ray &, const Eigen::MatrixXd &, const Eigen::MatrixXi &, std::vector<int> &); \
double approximate(const Ray &, const Triangle &); \
}
DECLARE_ADAPTER(baseline)
DECLARE_ADAPTER(shewchuk)
#undef DECLARE_ADAPTER
