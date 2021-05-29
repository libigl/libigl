function(igl_install module_name)
    if(NOT LIBIGL_INSTALL)
        return()
    endif()

    # Check if category is `copyleft` or `restricted`
    if(${module_name} MATCHES "^igl_copyleft")
        set(suffix "_copyleft")
    elseif(${module_name} MATCHES "^igl_restricted")
        set(suffix "_restricted")
    else()
        set(suffix "")
    endif()

    ########################
    # Install CMake target #
    ########################

    set_property(TARGET ${module_name} PROPERTY EXPORT_NAME ${module_export})

    include(GNUInstallDirs)
    install(TARGETS ${module_name}
        EXPORT LibiglTargets${suffix}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                COMPONENT LibiglRuntime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT          LibiglRuntime
                NAMELINK_COMPONENT LibiglDevelopment
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT LibiglRuntime
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
                COMPONENT LibiglDevelopment
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    ###################
    # Install headers #
    ###################

    # TODO: When moving module definition to their own CMakeLists.txt, we should
    # refactor this to use the folder where the target was defined (via the
    # target property SOURCE_DIR).
    set(target_include_dir ${PROJECT_SOURCE_DIR}/include)
    foreach(source_path IN ITEMS ${ARGN})
        # Filter out .cpp files in "static lib" mode
        if(LIBIGL_USE_STATIC_LIBRARY)
            get_filename_component(extension ${source_path} LAST_EXT)
            if(extension MATCHES ".cpp")
                continue()
            endif()
        endif()

        # Compute relative path to copy
        get_filename_component(source_directory ${source_path} DIRECTORY)
        file(RELATIVE_PATH source_subdir ${target_include_dir} ${source_directory})

        # Create install rule to copy specific header
        install(
            FILES ${source_path}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${source_subdir}
        )
    endforeach()
endfunction()
