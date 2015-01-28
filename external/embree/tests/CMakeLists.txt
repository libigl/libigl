## ======================================================================== ##
## Copyright 2009-2013 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_SOURCE_DIR}/../tutorials/common)

IF (__XEON__)

  IF (TARGET_SSE41)
    ADD_DEFINITIONS(-D__TARGET_SSE41__)
  ENDIF()

  IF (TARGET_AVX)
    ADD_DEFINITIONS(-D__TARGET_AVX__)
  ENDIF()

  IF (TARGET_AVX2)
    ADD_DEFINITIONS(-D__TARGET_AVX2__)
  ENDIF()

  ADD_EXECUTABLE(verify verify.cpp)
  TARGET_LINK_LIBRARIES(verify sys embree)

  ADD_EXECUTABLE(benchmark benchmark.cpp)
  TARGET_LINK_LIBRARIES(benchmark sys embree)

ELSE ()

  INCLUDE (icc_xeonphi)

  ADD_EXECUTABLE(verify_xeonphi verify.cpp)
  TARGET_LINK_LIBRARIES(verify_xeonphi sys_xeonphi embree_xeonphi.so)
  ADD_DEPENDENCIES(verify_xeonphi embree_xeonphi)

  ADD_EXECUTABLE(benchmark_xeonphi benchmark.cpp)
  TARGET_LINK_LIBRARIES(benchmark_xeonphi sys_xeonphi embree_xeonphi.so)
  ADD_DEPENDENCIES(benchmark_xeonphi embree_xeonphi)

ENDIF()
