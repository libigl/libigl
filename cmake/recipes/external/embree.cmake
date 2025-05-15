if(TARGET embree::embree)
    return()
endif()

message(STATUS "Third-party: creating target 'embree::embree'")

include(FetchContent)
FetchContent_Declare(
    embree
    GIT_REPOSITORY https://github.com/embree/embree.git
    GIT_TAG        v4.4.0
    GIT_SHALLOW    TRUE
)

# Set Embree's default options
option(EMBREE_ISPC_SUPPORT "Build Embree with support for ISPC applications." OFF)
option(EMBREE_TUTORIALS    "Enable to build Embree tutorials"                 OFF)
option(EMBREE_STATIC_LIB   "Build Embree as a static library."                ON)
set(EMBREE_TESTING_INTENSITY 0          CACHE STRING "Intensity of testing (0 = no testing, 1 = verify and tutorials, 2 = light testing, 3 = intensive testing.")
set(EMBREE_TASKING_SYSTEM    "INTERNAL" CACHE STRING "Selects tasking system")
set(EMBREE_MAX_ISA           "DEFAULT"  CACHE STRING "Selects highest ISA to support.")

# Ready to include embree's atrocious CMake
FetchContent_MakeAvailable(embree)

# Disable warnings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Embree's subgrid.h is known for causing array subscript out of bound
    # warning.  Embree dev claim the code is correct and it is a GCC bug
    # for misfiring warnings.  See https://github.com/embree/embree/issues/271
    #
    # The issue should be fixed for gcc 9.2.1 and later.
    target_compile_options(embree PRIVATE "-Wno-array-bounds")
endif()

# Now we need to do some juggling to propagate the include directory properties
# along with the `embree` target
include(GNUInstallDirs)
target_include_directories(embree SYSTEM INTERFACE
    "$<BUILD_INTERFACE:${embree_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>"
)

add_library(embree::embree ALIAS embree)

# Some order for IDEs
set_target_properties(embree PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(lexers PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(math PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(simd PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(sys PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(tasking PROPERTIES FOLDER "ThirdParty//Embree")
set_target_properties(uninstall PROPERTIES FOLDER "ThirdParty//Embree")
