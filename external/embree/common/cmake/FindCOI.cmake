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

FIND_PATH(COI_INCLUDE_PATH source/COIProcess_source.h /opt/intel/mic/coi/include /usr/include/intel-coi)
FIND_PATH(COI_HOST_LIBRARY_DIR libcoi_host.so /opt/intel/mic/coi/host-linux-release/lib/)
FIND_PATH(COI_DEV_LIBRARY_DIR libcoi_device.so /opt/mpss/*/sysroots/k1om-mpss-linux/usr/lib64/)

IF (COI_INCLUDE_PATH AND COI_HOST_LIBRARY_DIR AND COI_DEV_LIBRARY_DIR)
  SET(COI_FOUND TRUE)
  SET(COI_INCLUDE_PATHS ${COI_INCLUDE_PATH})
  SET(COI_HOST_LIBRARIES ${COI_HOST_LIBRARY_DIR}/libcoi_host.so)
  SET(COI_DEV_LIBRARIES ${COI_DEV_LIBRARY_DIR}/libcoi_device.so)
ELSE ()
  MESSAGE(FATAL_ERROR "COI installation not found")
ENDIF ()

MARK_AS_ADVANCED(COI_INCLUDE_PATH)
MARK_AS_ADVANCED(COI_HOST_LIBRARY_DIR)
MARK_AS_ADVANCED(COI_DEV_LIBRARY_DIR)

MARK_AS_ADVANCED(COI_INCLUDE_PATHS)
MARK_AS_ADVANCED(COI_HOST_LIBRARIES)
MARK_AS_ADVANCED(COI_DEV_LIBRARIES)
