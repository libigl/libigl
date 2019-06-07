################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option.
set(LIBIGL_EXTRA_OPTIONS TLS_VERIFY OFF)
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	list(APPEND LIBIGL_EXTRA_OPTIONS GIT_CONFIG advice.detachedHead=false)
endif()

# On CMake 3.6.3 and above, there is an option to use shallow clones of git repositories.
# The shallow clone option only works with real tags, not SHA1, so we use a separate option.
set(LIBIGL_BRANCH_OPTIONS)
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.3"))
	# Disabled for now until we can make sure that it has no adverse effects
	# (Downside is that the eigen mirror is huge again)
	# list(APPEND LIBIGL_BRANCH_OPTIONS GIT_SHALLOW 1)
endif()

option(LIBIGL_SKIP_DOWNLOAD "Skip downloading external libraries" OFF)

# Shortcut functions
function(igl_download_project_aux name source)
	if(NOT LIBIGL_SKIP_DOWNLOAD)
		download_project(
			PROJ         ${name}
			SOURCE_DIR   "${source}"
			DOWNLOAD_DIR "${LIBIGL_EXTERNAL}/.cache/${name}"
			QUIET
			${LIBIGL_EXTRA_OPTIONS}
			${ARGN}
		)
	endif()
endfunction()

function(igl_download_project name)
	igl_download_project_aux(${name} "${LIBIGL_EXTERNAL}/${name}" ${ARGN})
endfunction()

################################################################################

## CGAL
function(igl_download_cgal)
	igl_download_project(cgal
		GIT_REPOSITORY https://github.com/CGAL/cgal.git
		GIT_TAG        f7c3c8212b56c0d6dae63787efc99093f4383415
	)
endfunction()

## CoMISo
function(igl_download_comiso)
	igl_download_project(CoMISo
		GIT_REPOSITORY https://github.com/libigl/CoMISo.git
		GIT_TAG        1f9618cf9b7bd77370d817976470d59091928606
	)
endfunction()

## Cork
function(igl_download_cork)
	igl_download_project(cork
		GIT_REPOSITORY https://github.com/libigl/cork.git
		GIT_TAG        27ad8a285838f5a480d856429e39d3d56d4338f9
	)
endfunction()

## Eigen
set(LIBIGL_EIGEN_VERSION 3.2.10 CACHE STRING "Default version of Eigen used by libigl.")
function(igl_download_eigen)
	igl_download_project(eigen
		GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
		GIT_TAG        ${LIBIGL_EIGEN_VERSION}
		${LIBIGL_BRANCH_OPTIONS}
	)
endfunction()

## Embree
function(igl_download_embree)
	igl_download_project(embree
		GIT_REPOSITORY https://github.com/embree/embree.git
		GIT_TAG        v3.5.2
		${LIBIGL_BRANCH_OPTIONS}
	)
endfunction()

## glad
function(igl_download_glad)
	igl_download_project(glad
		GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
		GIT_TAG        09b4969c56779f7ddf8e6176ec1873184aec890f
	)
endfunction()

## GLFW
function(igl_download_glfw)
	igl_download_project(glfw
		GIT_REPOSITORY https://github.com/glfw/glfw.git
		GIT_TAG        3.3
		${LIBIGL_BRANCH_OPTIONS}
	)
endfunction()

## ImGui
function(igl_download_imgui)
	igl_download_project(imgui
		GIT_REPOSITORY https://github.com/ocornut/imgui.git
		GIT_TAG        v1.69
		${LIBIGL_BRANCH_OPTIONS}
	)
	igl_download_project(libigl-imgui
		GIT_REPOSITORY https://github.com/libigl/libigl-imgui.git
		GIT_TAG        07ecd3858acc71e70f0f9b2dea20a139bdddf8ae
	)
endfunction()

## pybind11
function(igl_download_pybind11)
	igl_download_project(pybind11
		GIT_REPOSITORY https://github.com/pybind/pybind11.git
		GIT_TAG        2d0507db43cd5a117f7843e053b17dffca114107
	)
endfunction()

## stb_image
function(igl_download_stb)
	igl_download_project(stb
		GIT_REPOSITORY https://github.com/libigl/libigl-stb.git
		GIT_TAG        cd0fa3fcd90325c83be4d697b00214e029f94ca3
	)
endfunction()

## TetGen
function(igl_download_tetgen)
	igl_download_project(tetgen
		GIT_REPOSITORY https://github.com/jdumas/tetgen.git
		GIT_TAG        c63e7a6434652b8a2065c835bd9d6d298db1a0bc
	)
endfunction()

## TinyXML
function(igl_download_tinyxml2)
	igl_download_project(tinyxml2
		GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
		GIT_TAG        d175e9de0be0d4db75d0a8cf065599a435a87eb6
	)
endfunction()

## Triangle
function(igl_download_triangle)
	igl_download_project(triangle
		GIT_REPOSITORY https://github.com/libigl/triangle.git
		GIT_TAG        d284c4a843efac043c310f5fa640b17cf7d96170
	)
endfunction()

## Catch2
function(igl_download_catch2)
	igl_download_project(catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        03d122a35c3f5c398c43095a87bc82ed44642516
	)
endfunction()

## Predicates
function(igl_download_predicates)
	igl_download_project(predicates
		GIT_REPOSITORY https://github.com/libigl/libigl-predicates.git
		GIT_TAG        4c57c1d3f31646b010d1d58bfbe201e75c2b2ad8
	)
endfunction()

################################################################################

## Test data
function(igl_download_test_data)
	igl_download_project_aux(test_data
		"${LIBIGL_EXTERNAL}/../tests/data"
		GIT_REPOSITORY https://github.com/libigl/libigl-tests-data
		GIT_TAG        adc66cabf712a0bd68ac182b4e7f8b5ba009c3dd
	)
endfunction()

## Tutorial data
function(igl_download_tutorial_data)
	igl_download_project_aux(tutorial_data
		"${LIBIGL_EXTERNAL}/../tutorial/data"
		GIT_REPOSITORY https://github.com/libigl/libigl-tutorial-data
		GIT_TAG        5c6a1ea809c043d71e5595775709c15325a7158c
	)
endfunction()

