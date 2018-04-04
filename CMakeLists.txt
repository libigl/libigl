#libigl2

add_library(libigl2 INTERFACE)

target_link_libraries(libigl2 INTERFACE eigen3 ethos)

target_include_directories(libigl2 INTERFACE
    include
)

target_compile_definitions(libigl2 INTERFACE
    _USE_MATH_DEFINES
)
