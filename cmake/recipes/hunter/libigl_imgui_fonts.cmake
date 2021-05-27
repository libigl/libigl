get_filename_component(current_filename ${CMAKE_CURRENT_LIST_FILE} NAME_WLE)
message(FATAL_ERROR
    "${current_filename} is not set up to be used with Hunter. "
    "Please disable the libigl modules depending on it, or set HUNTER_ENABLED to OFF."
)
