#
# Try to find a MATLAB installation
# Once done this will define
#
# MATLAB_FOUND
# MATLAB_INCLUDE_DIR
# MATLAB_LIBRARIES
#

FIND_PATH(MATLAB_INCLUDE_DIR engine.h
      PATHS
      /Applications/MATLAB_R2012b.app/extern/include
      /Applications/MATLAB_R2013a.app/extern/include
      /Applications/MATLAB_R2013b.app/extern/include
      /Applications/MATLAB_R2014a.app/extern/include
      /Applications/MATLAB_R2014b.app/extern/include
      NO_DEFAULT_PATH)

FIND_LIBRARY(MATLAB_LIBRARY1 eng
  PATHS
    /Applications/MATLAB_R2012b.app/bin/maci64
    /Applications/MATLAB_R2013a.app/bin/maci64
    /Applications/MATLAB_R2013b.app/bin/maci64
    /Applications/MATLAB_R2014a.app/bin/maci64
    /Applications/MATLAB_R2014b.app/bin/maci64
  PATH_SUFFIXES`
    a
    lib64
    lib
    dylib
    NO_DEFAULT_PATH
)

FIND_LIBRARY(MATLAB_LIBRARY2 mx
  PATHS
    /Applications/MATLAB_R2012b.app/bin/maci64/
    /Applications/MATLAB_R2013a.app/bin/maci64/
    /Applications/MATLAB_R2013b.app/bin/maci64/
    /Applications/MATLAB_R2014a.app/bin/maci64/
    /Applications/MATLAB_R2014b.app/bin/maci64/
  PATH_SUFFIXES
    a
    lib64
    lib
    dylib
    NO_DEFAULT_PATH
)

if(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARY1 AND MATLAB_LIBRARY2)
	message(STATUS "Found MATLAB: ${MATLAB_INCLUDE_DIR}")
  set(MATLAB_LIBRARIES ${MATLAB_LIBRARY1} ${MATLAB_LIBRARY2})
else(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARY1 AND MATLAB_LIBRARY2)
	message(FATAL_ERROR "could NOT find MATLAB")
endif(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARY1 AND MATLAB_LIBRARY2)
