set(TEST_SRC_FILES
	avg_edge_length.cpp
	bbw.cpp
	boundary_loop.cpp
	copyleft
	cotmatrix.cpp
	cotmatrix_entries.cpp
	cut_to_disk.cpp
	decimate.cpp
	dirname.cpp
	doublearea.cpp
	edge_flaps.cpp
	edge_lengths.cpp
	guess_extension.cpp
	is_edge_manifold.cpp
	is_symmetric.cpp
	ismember.cpp
	list_to_matrix.cpp
	main.cpp
	mosek
	pathinfo.cpp
	per_face_normals.cpp
	qslim.cpp
	readDMAT.cpp
	readOBJ.cpp
	readOFF.cpp
	seam_edges.cpp
	setdiff.cpp
	slice.cpp
	slice_into.cpp
	sort.cpp
	squared_edge_lengths.cpp
	tet_tet_adjacency.cpp
	triangle_triangle_adjacency.cpp
	unique.cpp
	upsample.cpp
)

add_executable(igl_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_tests igl::core gtest_main)
add_test(NAME run_igl_tests COMMAND igl_tests)
