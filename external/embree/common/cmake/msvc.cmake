## ======================================================================== ##
## Copyright 2009-2014 Intel Corporation                                    ##
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

# the /QxXXX flags are meant for the Intel Compiler, MSVC ignores them
SET(FLAGS_SSE2  "/QxSSE2")
SET(FLAGS_SSE3  "/QxSSE3")
SET(FLAGS_SSSE3 "/QxSSEE3")
SET(FLAGS_SSE41 "/DCONFIG_SSE41 /QxSSE4.1")
SET(FLAGS_SSE42 "/DCONFIG_SSE42 /QxSSE4.2")
SET(FLAGS_AVX   "/arch:AVX /DCONFIG_AVX")
# Intel Compiler 15, Update 1 unfortunately cannot handle /arch:AVX2
IF (COMPILER STREQUAL "ICC") # for scripts/regression.py to work with ICC
  SET(FLAGS_AVX2  "/DCONFIG_AVX2 /QxCORE-AVX2")
ELSE()
  SET(FLAGS_AVX2  "/arch:AVX2 /DCONFIG_AVX2 /QxCORE-AVX2")
ENDIF()
SET(FLAGS_AVX512 "")

SET(ADDITIONAL_CXX_FLAGS "/Ox /fp:fast /Qpar /Oi /Gy /GR- /MP")

SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} /MP")
SET(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS}")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${ADDITIONAL_CXX_FLAGS}")
# use static runtime library
string (REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
# remove define NDEBUG and instead set define DEBUG for config RelWithDebInfo
string (REPLACE "NDEBUG" "DEBUG" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
