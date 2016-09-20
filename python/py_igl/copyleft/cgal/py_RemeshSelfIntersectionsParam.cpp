py::class_<igl::copyleft::cgal::RemeshSelfIntersectionsParam > RemeshSelfIntersectionsParam(m, "RemeshSelfIntersectionsParam");

RemeshSelfIntersectionsParam
.def("__init__", [](igl::copyleft::cgal::RemeshSelfIntersectionsParam &m)
{
    new (&m) igl::copyleft::cgal::RemeshSelfIntersectionsParam();
    m.detect_only = false;
    m.first_only = false;
    m.stitch_all = false;
})
.def_readwrite("detect_only", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::detect_only)
.def_readwrite("first_only", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::first_only)
.def_readwrite("stitch_all", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::stitch_all)
;
