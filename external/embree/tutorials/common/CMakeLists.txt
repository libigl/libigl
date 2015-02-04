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

SET (__XEON__ 1)

# always add "image", because it has own logic to handle BUILD_TUTORIALS* flags 
ADD_SUBDIRECTORY(image)

IF (BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC)
  ADD_SUBDIRECTORY(lexers)
  ADD_SUBDIRECTORY(transport)
  ADD_SUBDIRECTORY(tutorial)
  
  IF (XEON_PHI_ISA)
    SET (__XEON__ 0)
    ADD_SUBDIRECTORY(tutorial tutorial_xeonphi)
    ADD_SUBDIRECTORY(transport_coi/host)
    ADD_SUBDIRECTORY(transport_coi/device)
  ENDIF()
ENDIF()
