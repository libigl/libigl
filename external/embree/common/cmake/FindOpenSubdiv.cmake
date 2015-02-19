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

FIND_PATH( OPENSUBDIV_INCLUDE_DIR NAMES opensubdiv/version.h )
FIND_LIBRARY( OPENSUBDIV_LIBRARY NAMES osdCPU)

SET(OPENSUBDIV_LIBRARIES ${OPENSUBDIV_LIBRARY} )
SET(OPENSUBDIV_INCLUDE_DIRS ${OPENSUBDIV_INCLUDE_DIR} )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENSUBDIV DEFAULT_MSG OPENSUBDIV_INCLUDE_DIR OPENSUBDIV_LIBRARY)

MARK_AS_ADVANCED(OPENSUBDIV_INCLUDE_DIR OPENSUBDIV_LIBRARY)
