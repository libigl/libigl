cmake_minimum_required(VERSION 3.1)
project(libigl_external)

include(ExternalProject)

set(CGAL_DIR "${CMAKE_CURRENT_LIST_DIR}/cgal")
ExternalProject_Add(CGAL_Project
	PREFIX "${CGAL_DIR}"

	GIT_REPOSITORY "https://github.com/CGAL/cgal.git"
	GIT_TAG "releases/CGAL-4.11"

	UPDATE_COMMAND ""
	PATCH_COMMAND ""
 
	CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF "-DCMAKE_INSTALL_PREFIX=${CGAL_DIR}/cgal/" -DCMAKE_INSTALL_LIBDIR=lib

	TEST_COMMAND ""
)
