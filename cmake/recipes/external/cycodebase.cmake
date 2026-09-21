if(TARGET cycodebase::cycodebase)
    return()
endif()

FetchContent_Declare(
  cyCodeBase
  GIT_REPOSITORY https://github.com/cemyuksel/cyCodeBase/
  GIT_TAG e36f3cffca65eb12a8a071f0443128b7de6ed75d
)
FetchContent_Populate(cyCodeBase)
add_library(cyCodeBase_interface INTERFACE)
target_include_directories(cyCodeBase_interface INTERFACE ${cycodebase_SOURCE_DIR})

if(NOT (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i[3-6]86"))
  target_compile_definitions(cyCodeBase_interface INTERFACE CY_NO_INTRIN_H)
endif()

add_library(cyCodeBase::cyCodeBase ALIAS cyCodeBase_interface)
