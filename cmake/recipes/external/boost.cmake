if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'...")

include(FetchContent)
FetchContent_Declare(
    boost-cmake
    GIT_REPOSITORY https://github.com/Orphis/boost-cmake.git
    GIT_TAG 70b12f62da331dd402b78102ec8f6a15d59a7af9
)

set(PREVIOUS_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(OLD_CMAKE_POSITION_INDEPENDENT_CODE ${CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This guy will download boost using FetchContent
FetchContent_MakeAvailable(boost-cmake)

set(CMAKE_POSITION_INDEPENDENT_CODE ${OLD_CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_CXX_FLAGS "${PREVIOUS_CMAKE_CXX_FLAGS}")

