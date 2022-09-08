if(TARGET Eigen3::Eigen)
    return()
endif()

message(STATUS "Third-party: creating target 'Eigen3::Eigen'")

include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG tags/3.3.7
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
endif()

add_library(Eigen3_Eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)

include(GNUInstallDirs)
target_include_directories(Eigen3_Eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

option(EIGEN_MPL2_ONLY "Compile Eigen with MPL2 code only" OFF)
if(EIGEN_MPL2_ONLY)
    target_compile_definitions(Eigen3_Eigen INTERFACE EIGEN_MPL2_ONLY)
endif()

# On Windows, enable natvis files to improve debugging experience
if(WIN32 AND eigen_SOURCE_DIR)
    target_sources(Eigen3_Eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}/debug/msvc/eigen.natvis>)
endif()

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME eigen)
set_target_properties(Eigen3_Eigen PROPERTIES EXPORT_NAME Eigen)
install(DIRECTORY ${eigen_SOURCE_DIR}/Eigen DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS Eigen3_Eigen EXPORT EigenTargets)
install(EXPORT EigenTargets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/eigen NAMESPACE Eigen3::)
