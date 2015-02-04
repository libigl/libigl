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


# additional parameters (beyond the name) are treated as additional dependencies
# if ADDITIONAL_LIBRARIES is set these will be included during linking

MACRO (ADD_TUTORIAL TUTORIAL_NAME)

IF (__XEON__)

IF (BUILD_TUTORIALS)
  ADD_EXECUTABLE(${TUTORIAL_NAME} ${TUTORIAL_NAME}.cpp ${TUTORIAL_NAME}_device.cpp ${ARGN})
  TARGET_LINK_LIBRARIES(${TUTORIAL_NAME} embree tutorial image transport tutorial_device ${ADDITIONAL_LIBRARIES})
  SET_PROPERTY(TARGET ${TUTORIAL_NAME} PROPERTY FOLDER tutorials/single)
ENDIF()

IF (BUILD_TUTORIALS_ISPC)
  ADD_ISPC_EXECUTABLE(${TUTORIAL_NAME}_ispc ${TUTORIAL_NAME}.cpp ${TUTORIAL_NAME}_device.ispc)
  TARGET_LINK_LIBRARIES(${TUTORIAL_NAME}_ispc embree tutorial image transport tutorial_device_ispc)
  SET_PROPERTY(TARGET ${TUTORIAL_NAME}_ispc PROPERTY FOLDER tutorials/ispc)
ENDIF()

ELSE (__XEON__)

IF (BUILD_TUTORIALS)
  IF (__HOST__)
    ADD_EXECUTABLE(${TUTORIAL_NAME}_xeonphi ${TUTORIAL_NAME}.cpp)
    TARGET_LINK_LIBRARIES(${TUTORIAL_NAME}_xeonphi tutorial image transport_host)
    SET_PROPERTY(TARGET ${TUTORIAL_NAME}_xeonphi PROPERTY FOLDER tutorials/xeonphi)
  ELSE()
    INCLUDE (icc_xeonphi)
    ADD_ISPC_EXECUTABLE(${TUTORIAL_NAME}_xeonphi_device ${TUTORIAL_NAME}_device.ispc)
    TARGET_LINK_LIBRARIES(${TUTORIAL_NAME}_xeonphi_device embree_xeonphi transport_device tutorial_xeonphi_device_ispc)
    SET_PROPERTY(TARGET ${TUTORIAL_NAME}_xeonphi_device PROPERTY FOLDER tutorials/xeonphi)
  ENDIF()
ENDIF()

ENDIF (__XEON__)


ENDMACRO ()