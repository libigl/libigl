cmake_minimum_required(VERSION 3.1)
project(libigl_external)

include(ExternalProject)

set(CGAL_DIR "${CMAKE_CURRENT_LIST_DIR}/cgal")
ExternalProject_Add(CGAL_Project
	PREFIX "${CGAL_DIR}"

	URL "https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.11/CGAL-4.11.zip"
	URL_MD5 320d70f7fd2a98537d20e21eeaea5a25

	UPDATE_COMMAND ""
	PATCH_COMMAND ""
 
	CMAKE_ARGS
    -DCMAKE_PREFIX_PATH=${CGAL_DIR}/auxiliary/gmp
    -DBUILD_SHARED_LIBS=OFF
    "-DCMAKE_INSTALL_PREFIX=${CGAL_DIR}/"
    -DCMAKE_INSTALL_LIBDIR=lib

	TEST_COMMAND ""
)

ExternalProject_Add_Step(
	CGAL_Project download_auxiliary
	COMMAND ${CMAKE_COMMAND} -DCGAL_DIR=${CGAL_DIR} -P ${CMAKE_CURRENT_LIST_DIR}/download.cmake
	DEPENDEES downdload
	DEPENDERS configure
)
