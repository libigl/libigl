if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'...")

cmake_minimum_required(VERSION 3.24) # Ensure modern FetchContent features
project(BoostFetchExample)

include(FetchContent)

# Define the Boost library to fetch
FetchContent_Declare(
    Boost
    URL https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
    URL_HASH MD5=ac857d73bb754b718a039830b07b9624
)
# Fetch Boost
FetchContent_MakeAvailable(Boost)


# Ensure Boost paths are set before CGAL
set(Boost_INCLUDE_DIR ${boost_SOURCE_DIR})
set(Boost_LIBRARY_DIR ${boost_BINARY_DIR})

# Add Boost libraries needed for your project
set(BOOST_LIBRARIES
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
    multiprecision
    program_options
    timer
    random
    serialization
    system
    thread
    type_erasure
  )

foreach(lib IN LISTS BOOST_LIBRARIES)
    add_library(boost_${lib} INTERFACE)
    target_include_directories(boost_${lib} INTERFACE ${Boost_SOURCE_DIR})
    target_link_libraries(boost_${lib} INTERFACE Boost::${lib})
endforeach()
