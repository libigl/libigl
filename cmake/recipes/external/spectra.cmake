if(TARGET Eigen3::Eigen)
  # Try to trick Spectra to avoid its find_package(Eigen3) getting confused.
  set(Eigen3_FOUND TRUE CACHE BOOL "" FORCE)
endif()
add_custom_target(eigen)


if(TARGET spectra::spectra)
  return()
endif()
include(FetchContent)

message(STATUS "Third-party: creating target 'spectra::spectra'")

FetchContent_Declare(
  Spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra/
  GIT_TAG v1.0.1
)
FetchContent_MakeAvailable(Spectra)

add_library(spectra::spectra ALIAS Spectra)
