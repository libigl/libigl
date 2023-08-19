function(igl_add_tutorial name)
    # Only create tutorial if dependencies have been enabled in the current build
    foreach(module IN ITEMS ${ARGN})
        if(NOT TARGET ${module})
            return()
        endif()
    endforeach()

    message(STATUS "Creating libigl tutorial: ${name}")
    # get all cpp files in ${CMAKE_CURRENT_SOURCE_DIR}/${name}/
    file(GLOB SRCFILES ${CMAKE_CURRENT_SOURCE_DIR}/${name}/*.cpp)
    add_executable(${name} ${SRCFILES})
    target_link_libraries(${name} PRIVATE
        igl::core
        igl::tutorial_data
        ${ARGN}
    )

    set_target_properties(${name} PROPERTIES FOLDER Libigl_Tutorials)
endfunction()
