function(igl_add_tutorial name)
    add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}/main.cpp)
    target_link_libraries(${name} PRIVATE
        igl::core
        igl::tutorial_data
        ${ARGN}
    )
endfunction()
