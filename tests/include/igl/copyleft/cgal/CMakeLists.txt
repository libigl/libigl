set(TEST_SRC_FILES
	CSGTree.cpp
	hausdorff.cpp
	main.cpp
	order_facets_around_edges.cpp
	outer_facet.cpp
	outer_hull.cpp
	peel_outer_hull_layers.cpp
	points_inside_component.cpp
	remesh_self_intersections.cpp
)

add_executable(igl_cgal_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_cgal_tests igl::core igl::cgal gtest_main)
add_test(NAME run_igl_cgal_tests COMMAND igl_cgal_tests)
