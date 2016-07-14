py::class_<RotationList>(m, "RotationList")
    .def(py::init<>())
    .def(py::init<size_t>())
    .def("pop_back", &RotationList::pop_back)
    /* There are multiple versions of push_back(), etc. Select the right ones. */
    .def("append", (void (RotationList::*)(const Eigen::Quaterniond &)) &RotationList::push_back)
    .def("back", (Eigen::Quaterniond &(RotationList::*)()) &RotationList::back)
    .def("__len__", [](const RotationList &v) { return v.size(); })
    .def("__getitem__", [](const RotationList &v, int b) { return v.at(b); })
    .def("__setitem__", [](RotationList &v, int b, Eigen::Quaterniond &c) { return v.at(b) = c; })
    .def("__iter__", [](RotationList &v) {
       return py::make_iterator(v.begin(), v.end());
}, py::keep_alive<0, 1>());
