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

SET(ISPC_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${ISPC_TARGET_DIR})

SET(ISPC_DIR $ENV{ISPC_DIR} CACHE STRING "Path to ISPC root dir.")
MARK_AS_ADVANCED(ISPC_DIR)

IF (${TARGET_SPMD} STREQUAL "knc")
  INCLUDE(build_icc_mic)
  SET(IVL_TARGET "knc")
  SET (ISPC_TARGET generic-16 --emit-c++ --c++-include-file=${ISPC_DIR}/examples/intrinsics/knc.h)
  INCLUDE_DIRECTORIES(/opt/intel/mic/coi/include)
  SET (ISPC_TARGET_EXT "cpp")
  SET (ISPC_TARGET_DEFINE __MIC__)
  ADD_DEFINITIONS(-D__ISPC_TARGET_MIC__)
ELSEIF (${TARGET_SPMD} STREQUAL "avx")
  SET (ISPC_TARGET avx)
  SET (ISPC_TARGET_EXT o)
  SET (ISPC_TARGET_DEFINE __AVX__)
  ADD_DEFINITIONS(-D__ISPC_TARGET_AVX__)
ELSEIF (${TARGET_SPMD} STREQUAL "sse")
  SET (ISPC_TARGET sse4)
  SET (ISPC_TARGET_EXT o)
  SET (ISPC_TARGET_DEFINE __SSE4__)
  ADD_DEFINITIONS(-D__ISPC_TARGET_SSE__)
ELSE ()
  MESSAGE(FATAL_ERROR "Unknown backend specified: " ${TARGET_SPMD})
ENDIF ()

#######################################################
# rule to compile .ispc file to .dev.cpp and .dev.h
# Note that this is only a *rule*, so as long as no other
# library or executable actually *uses* the respective
# .dev.cpp file nothing at all will happen.
#######################################################
MACRO(ISPC_COMPILE _filename _dest_base)

  GET_FILENAME_COMPONENT(_basename ${_filename} NAME_WE)
  GET_FILENAME_COMPONENT(_dirname ${_filename} PATH)
  
  IF("${_dirname}" STREQUAL "")
    SET(_full_dir_name ${_dest_base})
  ELSE("${_dirname}" STREQUAL "")
    SET(_full_dir_name ${_dest_base}/${_dirname})
  ENDIF("${_dirname}" STREQUAL "")
  
  SET(_out_base ${_full_dir_name}/${_basename})
  SET(_out_base2 ${ISPC_TARGET_DIR}/${_basename})

  IF (EXISTS ${_out_base}.dev.idep)
    FILE(READ ${_out_base}.dev.idep contents)
    STRING(REGEX REPLACE " " ";" contents "${contents}")
    STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
    STRING(REGEX REPLACE "\n" ";" contents "${contents}")
    SET(deps "")
    FOREACH(dep ${contents})
      IF (EXISTS ${dep})
        SET(deps ${deps} ${dep})
      ENDIF (EXISTS ${dep})
    ENDFOREACH(dep ${contents})
  ELSE (EXISTS ${_out_base}.dev.idep )
    SET(deps "")
  ENDIF ()

  ADD_CUSTOM_COMMAND(
      OUTPUT ${_out_base}.dev.${ISPC_TARGET_EXT} ${_out_base2}_ispc.h
      COMMAND mkdir -p ${_full_dir_name} \; ispc  
      -I ${CMAKE_CURRENT_SOURCE_DIR} 
      -D${ISPC_TARGET_DEFINE}
      --pic
      # FIXME                 -g      
      -O1 
      --target ${ISPC_TARGET}
      --wno-perf
      --opt=fast-math
      --opt=force-aligned-memory
      -h ${_out_base2}_ispc.h
      -o ${_out_base}.dev.${ISPC_TARGET_EXT}
      -MMM  ${_out_base}.dev.idep 
      ${CMAKE_CURRENT_SOURCE_DIR}/${_filename} 
      \;
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_filename}
      ${deps})

ENDMACRO(ISPC_COMPILE)

FOREACH(ispc_src ${ISPC_SOURCES})
  GET_FILENAME_COMPONENT(_basename ${ispc_src} NAME_WE)
  GET_FILENAME_COMPONENT(_dirname ${ispc_src} PATH)
  ISPC_COMPILE(${ispc_src} ${ISPC_TARGET_DIR})
  SET(SOURCES ${SOURCES} "${ISPC_TARGET_DIR}/${_dirname}/${_basename}.dev.${ISPC_TARGET_EXT}")
ENDFOREACH()

