m.def("remesh_self_intersections", []
(
	const Eigen::MatrixXd& V,	
	const Eigen::MatrixXi& F,
	const igl::copyleft::cgal::RemeshSelfIntersectionsParam& params,
	Eigen::MatrixXd& VV,
	Eigen::MatrixXi& FF,
	Eigen::MatrixXi& IF,
	Eigen::MatrixXi& J,
	Eigen::MatrixXi& IM
)
{
	assert_is_VectorX("J", J);
	assert_is_VectorX("IM", IM);
	Eigen::VectorXi Jt;
	Eigen::VectorXi IMt;
	igl::copyleft::cgal::remesh_self_intersections(V, F, params, VV, FF, IF, Jt, IMt);
	J = Jt;
	IM = IMt;
}, __doc_igl_copyleft_cgal_remesh_self_intersections,
py::arg("V"), py::arg("F"), py::arg("params"), py::arg("VV")
, py::arg("FF"), py::arg("IF"), py::arg("J"), py::arg("IM")); 