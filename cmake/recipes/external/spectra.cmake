if(TARGET spectra::spectra)
  return()
endif()
include(FetchContent)

message(STATUS "Third-party: creating target 'spectra::spectra'")

# Use fork because yixuan/spectra struggles to find Eigen3
FetchContent_Declare(
  Spectra
  GIT_REPOSITORY https://github.com/alecjacobson/spectra/
  GIT_TAG bbdc521b70a733c52ebfc0ac1484c82e13c3d140
)
FetchContent_MakeAvailable(Spectra)

add_library(spectra::spectra ALIAS Spectra)
