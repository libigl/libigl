target_link_libraries(libigl2 eigen3 ethos)

target_include_directories(libigl2 PUBLIC
    include
)

target_compile_definitions(libigl2 PUBLIC
    _USE_MATH_DEFINES
)
