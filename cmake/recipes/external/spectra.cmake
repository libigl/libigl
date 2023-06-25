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
