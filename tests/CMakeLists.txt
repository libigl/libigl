# This file is part of libigl, a simple c++ geometry processing library.
# 
# Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public License 
# v. 2.0. If a copy of the MPL was not distributed with this file, You can 
# obtain one at http://mozilla.org/MPL/2.0/.
#
# This file is based on PyMesh's unit test setup.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Set compiler
IF(APPLE)
    SET(CMAKE_C_COMPILER "clang")
    SET(CMAKE_CXX_COMPILER "clang++")
ELSEIF(UNIX)
    SET(CMAKE_C_COMPILER "gcc")
    SET(CMAKE_CXX_COMPILER "g++")
ELSEIF(WIN32)
    # TODO: not how the code would work on Windows.
ENDIF()

PROJECT(libigl_unit_tests)

INCLUDE(Settings.cmake)

# Process code in each subdirectories
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/include/igl)
