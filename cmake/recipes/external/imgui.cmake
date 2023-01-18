if(TARGET imgui::imgui)
    return()
endif()

message(STATUS "Third-party: creating target 'imgui::imgui'")

include(FetchContent)
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.85
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(imgui)

set(IMGUI_SRC
    "${imgui_SOURCE_DIR}/imgui.h"
    "${imgui_SOURCE_DIR}/imconfig.h"
    "${imgui_SOURCE_DIR}/imgui_internal.h"
    "${imgui_SOURCE_DIR}/imgui.cpp"
    "${imgui_SOURCE_DIR}/imgui_demo.cpp"
    "${imgui_SOURCE_DIR}/imgui_draw.cpp"
    "${imgui_SOURCE_DIR}/imgui_widgets.cpp"
    "${imgui_SOURCE_DIR}/imgui_tables.cpp"
    "${imgui_SOURCE_DIR}/imstb_rectpack.h"
    "${imgui_SOURCE_DIR}/imstb_textedit.h"
    "${imgui_SOURCE_DIR}/imstb_truetype.h"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.h"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.h"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp"
    "${imgui_SOURCE_DIR}/misc/cpp/imgui_stdlib.h"
    "${imgui_SOURCE_DIR}/misc/cpp/imgui_stdlib.cpp"
)

add_library(imgui ${IMGUI_SRC})
add_library(imgui::imgui ALIAS imgui)

# Include headers
target_include_directories(imgui PUBLIC "${imgui_SOURCE_DIR}")

# Compile definitions
target_compile_definitions(imgui PUBLIC
    IMGUI_IMPL_OPENGL_LOADER_GLAD
    IMGUI_DISABLE_OBSOLETE_FUNCTIONS # to check for obsolete functions
)

# Dependencies
include(glfw)
include(glad)
target_link_libraries(imgui PUBLIC glfw::glfw glad::glad)

set_target_properties(imgui PROPERTIES FOLDER ThirdParty)
