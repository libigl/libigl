if(TARGET CoMISo::CoMISo)
    return()
endif()

message(STATUS "Third-party: creating target 'CoMISo::CoMISo'")

include(FetchContent)
FetchContent_Declare(
    comiso
    GIT_REPOSITORY https://github.com/libigl/CoMISo.git
    GIT_TAG 536440e714f412e7ef6c0b96b90ba37b1531bb39
)

include(eigen)

FetchContent_MakeAvailable(comiso)

add_library(CoMISo::CoMISo ALIAS CoMISo)

# Copy .hh headers into a subfolder `CoMISo/`
file(GLOB_RECURSE INC_FILES "${comiso_SOURCE_DIR}/*.hh" "${comiso_SOURCE_DIR}/*.cc")
set(output_folder "${CMAKE_CURRENT_BINARY_DIR}/CoMISo/include/CoMISo")
message(VERBOSE "Copying CoMISo headers to '${output_folder}'")
foreach(filepath IN ITEMS ${INC_FILES})
    file(RELATIVE_PATH filename "${comiso_SOURCE_DIR}" ${filepath})
    configure_file(${filepath} "${output_folder}/${filename}" COPYONLY)
endforeach()

target_include_directories(CoMISo PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/CoMISo/include)

set_target_properties(CoMISo PROPERTIES FOLDER ThirdParty)
