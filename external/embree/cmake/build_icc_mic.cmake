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

#SET(COMPILER_TARGET "MIC" CACHE STRING "SSE version to use (SSSE3,SSE4.1,SSE4.2,AVX)" FORCE)

SET(__SSE__ 0)
SET(__AVX__ 0)

SET (CMAKE_CXX_COMPILER icpc)
SET (CMAKE_CXX_FLAGS "-mmic -restrict -Wall -fasm-blocks -fPIC")
SET (CMAKE_CXX_FLAGS_NOOPT "-O0 -DDEBUG")
SET (CMAKE_CXX_FLAGS_DEBUG "-g -w1 -O2 -DDEBUG ")

# -ansi-alias # ISPC not working with ansi-alias flag

SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mCG_lrb_num_threads=4 -mCG_lrb_num_threads=4 -fp-model fast -fimf-precision=low -fasm-blocks -no-inline-max-total-size -inline-factor=200 -fPIC  -fma  -restrict -no-prec-div -no-prec-sqrt  -mGLOB_default_function_attrs=\"use_vec_for_imul=on;use_fast_math=on;gather_scatter_loop_jknzd=on;gather_scatter_loop_unroll=2;use_gather_scatter_hint=on;c_lrb_avoid_vector_insts_with_rip=on;c_avoid_bank_conflicts=on;c_sch_nop_insertion=on;c_avoid_movz_and_movs=off;c_avoid_integer_ciscization=on;avoid_stall_between_all_insts=on;avoid_long_vector_ints=on;avoid_loads_with_extend=on;smart_mem_conflicts=on\"  -mP2OPT_hlo_prefetch=F ")

SET(CMAKE_LINKER icpc)
SET(CMAKE_CXX_LINK_EXECUTABLE	"<CMAKE_CXX_COMPILER> -mmic -static-intel -L/opt/intel/composerxe_mic/compiler/lib/mic -rdynamic -fPIC -L/opt/intel/mic/coi/device-linux-release/lib <LINK_FLAGS> -o <TARGET> ${CMAKE_GNULD_IMAGE_VERSION} <OBJECTS> <LINK_LIBRARIES> -zmuldefs")
SET(CMAKE_CXX_CREATE_SHARED_LIBRARY "<CMAKE_CXX_COMPILER> -mmic -static-intel  -L/opt/intel/composerxe_mic/compiler/lib/mic/ -L/opt/intel/mic/coi/device-linux-release  <LANGUAGE_COMPILE_FLAGS> <CMAKE_SHARED_MODULE_CXX_FLAGS> <LINK_FLAGS> <CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS> -o <TARGET> ${CMAKE_GNULD_IMAGE_VERSION} <OBJECTS> <LINK_LIBRARIES>   -zmuldefs")

SET(EXT "_knc")
