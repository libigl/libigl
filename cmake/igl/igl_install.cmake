function(igl_install module_name)
  if (NOT LIBIGL_INSTALL)
    return()
  endif ()

  # Check module name
  if (NOT ${module_name} MATCHES "^igl_")
    message(FATAL_ERROR "Libigl module name should start with 'igl_'")
  endif ()

  # extract suffix & component from module name
  if (${module_name} MATCHES "^igl_copyleft")
    set(suffix "-copyleft")
  elseif (${module_name} MATCHES "^igl_restricted")
    set(suffix "-restricted")
  else ()
    set(suffix "")
  endif ()

  string(REGEX REPLACE "^.*_?.+_" "-" component ${module_name})

  ########################
  # Install CMake target #
  ########################

  set(exports_name ${PROJECT_NAME}${suffix}${component}-targets)

  include(GNUInstallDirs)
  install(TARGETS ${module_name}
    EXPORT ${exports_name}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT LibiglRuntime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT LibiglRuntime
    NAMELINK_COMPONENT LibiglDevelopment
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT LibiglRuntime
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT LibiglDevelopment
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

  ################################
  # Install Config Package Files #
  ################################

  include(GNUInstallDirs)
  set(module_config_in "${PROJECT_SOURCE_DIR}/cmake/igl/libigl${suffix}${component}-config.cmake.in")
  set(module_config_out "${CMAKE_CURRENT_BINARY_DIR}/libigl${suffix}${component}-config.cmake")
  set(export_dest_dir "${CMAKE_INSTALL_LIBDIR}/cmake/libigl")

  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    "${module_config_in}"
    "${module_config_out}"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/libigl/cmake
  )
  install(FILES "${module_config_out}" DESTINATION "${export_dest_dir}")

  string(REPLACE "-" "_" namespace_suffix "${suffix}")
  install(EXPORT ${exports_name}
    DESTINATION ${export_dest_dir}
    NAMESPACE igl${namespace_suffix}::
    COMPONENT LibiglDevelopment)

  ###################
  # Install headers #
  ###################

  # TODO: When moving module definition to their own CMakeLists.txt, we should
  # refactor this to use the folder where the target was defined (via the
  # target property SOURCE_DIR).
  set(target_include_dir ${libigl_SOURCE_DIR}/include)
  foreach (source_path IN ITEMS ${ARGN})
    # Filter out .cpp files in "static lib" mode
    if (LIBIGL_USE_STATIC_LIBRARY)
      get_filename_component(extension ${source_path} LAST_EXT)
      if (extension MATCHES ".cpp")
        continue()
      endif ()
    endif ()

    # Compute relative path to copy
    get_filename_component(source_directory ${source_path} DIRECTORY)
    file(RELATIVE_PATH source_subdir ${target_include_dir} ${source_directory})

    # Create install rule to copy specific header
    install(
      FILES ${source_path}
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${source_subdir}
    )
  endforeach ()
endfunction()
