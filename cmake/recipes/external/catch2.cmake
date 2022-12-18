if (TARGET Catch2::Catch2)
  return()
endif ()

message(STATUS "Third-party: creating target 'Catch2::Catch2'")

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.8
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

include("${Catch2_SOURCE_DIR}/contrib/Catch.cmake")
