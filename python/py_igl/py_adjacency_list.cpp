m.def("adjacency_list", [](const Eigen::MatrixXi& F, std::vector<std::vector<int>>& A, bool sorted) {
    igl::adjacency_list(F, A, sorted);
}, py::arg("F"), py::arg("A"), py::arg("sorted")=false);
