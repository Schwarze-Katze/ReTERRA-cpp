cmake_minimum_required(VERSION 3.1)
project(lkhExec)
include_directories(${PROJECT_SOURCE_DIR}/INCLUDE)
add_definitions(-DTWO_LEVEL_TREE)
add_compile_options(-O3 -Wall -g -fcommon -flto)
file(GLOB LKHSRC "${PROJECT_SOURCE_DIR}/*.c")
add_executable(${PROJECT_NAME} ${LKHSRC})
target_link_libraries(${PROJECT_NAME} PRIVATE m)
