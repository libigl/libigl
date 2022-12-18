add_test(
  NAME
  COMMAND
  "${CMAKE_COMMAND}"
  --install "${CMAKE_CURRENT_BINARY_DIR}/${project_name}"
  --prefix "${test_install_dir}"
  --config $<CONFIG>
  --verbose)
