function(igl_add_tutorial name)
    # Only create tutorial if dependencies have been enabled in the current build
    foreach(module IN ITEMS ${ARGN})
        if(NOT TARGET ${module})
            return()
        endif()
    endforeach()

    message(STATUS "Creating libigl tutorial: ${name}")
    add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}/main.cpp)
    target_link_libraries(${name} PRIVATE
        igl::core
        igl::tutorial_data
        ${ARGN}
    )

    set_target_properties(${name} PROPERTIES FOLDER Libigl_Tutorials)
    if(LIBIGL_DELETE_TUTORIALS_IMMEDIATELY)
      message(STATUS "Adding command to remove ${name} immediately after build")
      add_custom_command(
       TARGET ${name}
       POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E remove $<TARGET_FILE:${name}>
       COMMENT "Removing ${name}.")
    endif()
endfunction()
