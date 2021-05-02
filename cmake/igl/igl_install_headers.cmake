function(igl_install_headers module_name)
    # TODO:
    # 1. Iterate over ARGN
    # 2. Get path relative to module_name's folder
    # 3. Filter out .cpp if using STATIC mode (temporary)
    # 4. Copy to destination using install() command
    #
    # install(
    #     FILES ${files_to_install}
    #     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/igl${subpath}
    # )
endfunction()
