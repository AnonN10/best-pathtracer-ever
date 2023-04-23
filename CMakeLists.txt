cmake_minimum_required(VERSION 3.24)

# set the project name
project(HelloWorld)

# add the executable
add_executable(HelloWorld main.cpp)

add_subdirectory(contrib/SDL)
add_subdirectory(contrib/bvh)
add_subdirectory(contrib/glm)

# libraries
target_link_libraries(${PROJECT_NAME} PRIVATE SDL2::SDL2 SDL2::SDL2main glm)

# Windows-only configurations
if (WIN32)
# copy the .dll file to the same folder as the executable
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:SDL2::SDL2>
    $<TARGET_FILE_DIR:${PROJECT_NAME}>
    VERBATIM)
endif()

# copy resources to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/bluenoise256.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR})