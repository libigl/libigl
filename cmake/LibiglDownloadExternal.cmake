################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	set(LIBIGL_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
	set(LIBIGL_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(igl_download_project name)
	download_project(
		PROJ         ${name}
		SOURCE_DIR   ${LIBIGL_EXTERNAL}/${name}
		DOWNLOAD_DIR ${LIBIGL_EXTERNAL}/.cache/${name}
		QUIET
		${LIBIGL_EXTRA_OPTIONS}
		${ARGN}
	)
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
		GIT_TAG        fea3ee0ba7d42ee3eca202d484e4fad855e4d6aa
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
function(igl_download_eigen)
	igl_download_project(eigen
		URL           http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
		URL_MD5       8ad10ac703a78143a4062c9bda9d8fd3
	)
endfunction()

## Embree
function(igl_download_embree)
	igl_download_project(embree
		URL            https://github.com/embree/embree/archive/v2.17.4.tar.gz
		URL_MD5        2038f3216b1d626e87453aee72c470e5
		# GIT_REPOSITORY https://github.com/embree/embree.git
		# GIT_TAG        cb61322db3bb7082caed21913ad14869b436fe78
	)
endfunction()

## glad
function(igl_download_glad)
	igl_download_project(glad
		GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
		GIT_TAG        71e35fe685a0cc160068a2f2f971c40b82d14af0
	)
endfunction()

## GLFW
function(igl_download_glfw)
	igl_download_project(glfw
		GIT_REPOSITORY https://github.com/glfw/glfw.git
		GIT_TAG        58cc4b2c5c2c9a245e09451437dd6f5af4d60c84
	)
endfunction()

## ImGui
function(igl_download_imgui)
	igl_download_project(imgui
		GIT_REPOSITORY https://github.com/ocornut/imgui.git
		GIT_TAG        bc6ac8b2aee0614debd940e45bc9cd0d9b355c86
	)
	igl_download_project(libigl-imgui
		GIT_REPOSITORY https://github.com/libigl/libigl-imgui.git
		GIT_TAG        a37e6e59e72fb07bd787dc7e90f72b9e1928dae7
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
		GIT_TAG        e671ceb3def5e7029a23de14c55dc16301ad4dab
	)
endfunction()

## TetGen
function(igl_download_tetgen)
	igl_download_project(tetgen
		GIT_REPOSITORY https://github.com/libigl/tetgen.git
		GIT_TAG        d2dcc33cb8551f16e302c8130ce12fa52992cd09
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
		GIT_TAG        d6761dd691e2e1318c83bf7773fea88d9437464a
	)
endfunction()

## Google test
function(igl_download_googletest)
	igl_download_project(googletest
		GIT_REPOSITORY https://github.com/google/googletest
		GIT_TAG        release-1.8.1
	)
endfunction()

################################################################################

## Test data
function(igl_download_test_data)
	set(IGL_TEST_DATA ${LIBIGL_EXTERNAL}/../tests/data)

	download_project(
		PROJ         test_data
		SOURCE_DIR   ${IGL_TEST_DATA}
		DOWNLOAD_DIR ${LIBIGL_EXTERNAL}/.cache/test_data
		QUIET
		GIT_REPOSITORY https://github.com/libigl/libigl-tests-data
		GIT_TAG        c81bb3b3db4cfd78bac6d359d845c45bc1059c9a
		${LIBIGL_EXTRA_OPTIONS}
	)
endfunction()

## Tutorial data
function(igl_download_tutorial_data)
	set(IGL_TUTORIAL_DATA ${LIBIGL_EXTERNAL}/../tutorial/data)

	download_project(
		PROJ         tutorial_data
		SOURCE_DIR   ${IGL_TUTORIAL_DATA}
		DOWNLOAD_DIR ${LIBIGL_EXTERNAL}/.cache/tutorial_data
		QUIET
		GIT_REPOSITORY https://github.com/libigl/libigl-tutorial-data
		GIT_TAG        5c6a1ea809c043d71e5595775709c15325a7158c
		${LIBIGL_EXTRA_OPTIONS}
	)
endfunction()

