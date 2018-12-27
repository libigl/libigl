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
function(igl_download_eigen)
	igl_download_project(eigen
		URL           http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
		URL_MD5       8ad10ac703a78143a4062c9bda9d8fd3
	)
endfunction()

## Embree
function(igl_download_embree)
	igl_download_project(embree
		URL            https://github.com/embree/embree/archive/v3.2.3.tar.gz
		URL_MD5        1868cda1c97d83d7a0b67b0b64b18cef
		# GIT_REPOSITORY https://github.com/embree/embree.git
		# GIT_TAG        cb61322db3bb7082caed21913ad14869b436fe78
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
		GIT_TAG        53c8c72c676ca97c10aedfe3d0eb4271c5b23dba
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
		GIT_TAG        edfa26e389060c21b9dd7812a0b19c00208b7224
	)
endfunction()

## TetGen
function(igl_download_tetgen)
	igl_download_project(tetgen
		GIT_REPOSITORY https://github.com/jdumas/tetgen.git
		GIT_TAG        63b4bdc5b947f9db75f01e0da36af54074ace5c9
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
		GIT_REPOSITORY https://github.com/jdumas/triangle.git
		GIT_TAG        2cd0672ff1f67f9f6bb8e556e84901293e637b76
	)
endfunction()

## Catch2
function(igl_download_catch2)
	igl_download_project(catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        03d122a35c3f5c398c43095a87bc82ed44642516
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
		GIT_TAG        adc66cabf712a0bd68ac182b4e7f8b5ba009c3dd
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

