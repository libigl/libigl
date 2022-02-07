if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'...")

include(FetchContent)
FetchContent_Declare(
    boost-cmake
    GIT_REPOSITORY https://github.com/Orphis/boost-cmake.git
    GIT_TAG 7f97a08b64bd5d2e53e932ddf80c40544cf45edf
)

set(PREVIOUS_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(OLD_CMAKE_POSITION_INDEPENDENT_CODE ${CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This guy will download boost using FetchContent
FetchContent_GetProperties(boost-cmake)
if(NOT boost-cmake_POPULATED)
    FetchContent_Populate(boost-cmake)
    # File lcid.cpp from Boost_locale.cpp doesn't compile on MSVC, so we exclude them from the default
    # targets being built by the project (only targets explicitly used by other targets will be built).
    add_subdirectory(${boost-cmake_SOURCE_DIR} ${boost-cmake_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ${OLD_CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_CXX_FLAGS "${PREVIOUS_CMAKE_CXX_FLAGS}")

# Set VS target folders
set(boost_modules
    container
    regex
    atomic
    exception
    chrono
    wave
    context
    coroutine
    date_time
    fiber
    filesystem
    graph
    iostreams
    locale
    log
    log_setup
    unit_test_framework
    math
    program_options
    timer
    random
    serialization
    system
    thread
    type_erasure
)
foreach(module IN ITEMS ${boost_modules})
    if(TARGET Boost_${module})
        set_target_properties(Boost_${module} PROPERTIES FOLDER ThirdParty/Boost)
    endif()
endforeach()

