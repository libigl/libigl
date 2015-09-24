# This file is part of libigl, a simple c++ geometry processing library.
# 
# Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public License 
# v. 2.0. If a copy of the MPL was not distributed with this file, You can 
# obtain one at http://mozilla.org/MPL/2.0/.
#
# This file is based on PyMesh's unit test setup.

# Include directories to search for source.
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
GET_FILENAME_COMPONENT(LIBIGL_PATH ${PROJECT_SOURCE_DIR} DIRECTORY)
INCLUDE_DIRECTORIES(${LIBIGL_PATH}/include/)

# Set build type.
SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)

# Create 64 bits binary.  32 bits support is dropped.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Os")
SET(CMAKE_LIBRARY_PATH /opt/local/lib ${CMAKE_LIBRARY_PATH})

# Set output directories
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
MAKE_DIRECTORY(${EXECUTABLE_OUTPUT_PATH})

LINK_DIRECTORIES(/opt/local/lib)

SET(CMAKE_MACOSX_RPATH ON)

# Set TEST_DIR definition
ADD_DEFINITIONS(-DTEST_DIR="${PROJECT_SOURCE_DIR}")

# Include current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Include Eigen
#SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(CMAKE_MODULE_PATH ${LIBIGL_PATH}/tutorial/cmake)
FIND_PACKAGE(Eigen REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS}/unsupported)
ADD_DEFINITIONS(-DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)

# Add googletest googlemock support
ADD_SUBDIRECTORY(${LIBIGL_PATH}/external/googletest/googlemock
    ${CMAKE_CURRENT_BINARY_DIR}/gtest)
SET(GTEST_BOTH_LIBRARIES gtest gmock)
INCLUDE_DIRECTORIES(${gmock_SOURCE_DIR})
INCLUDE_DIRECTORIES(${gmock_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR})
INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include)
