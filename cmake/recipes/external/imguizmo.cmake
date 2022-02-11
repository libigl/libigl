if(TARGET imguizmo::imguizmo)
    return()
endif()

message(STATUS "Third-party: creating target 'imguizmo::imguizmo'")

include(FetchContent)
FetchContent_Declare(
    imguizmo
    GIT_REPOSITORY https://github.com/CedricGuillemet/ImGuizmo.git
    GIT_TAG a23567269f6617342bcc112394bdad937b54b2d7
)
FetchContent_MakeAvailable(imguizmo)

set(IMGUIZMO_SRC
    "${imguizmo_SOURCE_DIR}/ImGuizmo.h"
    "${imguizmo_SOURCE_DIR}/ImGuizmo.cpp"
)

add_library(imguizmo ${IMGUIZMO_SRC})
add_library(imguizmo::imguizmo ALIAS imguizmo)

target_compile_features(imguizmo PUBLIC cxx_std_11)

target_include_directories(imguizmo PUBLIC "${imguizmo_SOURCE_DIR}")

include(imgui)
target_link_libraries(imguizmo PUBLIC imgui::imgui)

set_target_properties(imguizmo PROPERTIES FOLDER ThirdParty)
