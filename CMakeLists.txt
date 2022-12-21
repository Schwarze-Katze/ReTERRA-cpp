cmake_minimum_required(VERSION 3.1)
project(reTERRA-cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

#This project requires CGAL
find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()
include_directories(${CGAL_INCLUDE_DIRS})
link_libraries(${CGAL_LIBRARIES})

#CGAL requires Boost
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()


include_directories(
    ${PROJECT_SOURCE_DIR}/include

)
add_executable(ReTERRA
    ${PROJECT_SOURCE_DIR}/src/TERRAExec.cpp
    ${PROJECT_SOURCE_DIR}/src/TERRA.cpp
    ${PROJECT_SOURCE_DIR}/src/Voronoi.cpp
    ${PROJECT_SOURCE_DIR}/src/SetCover.cpp
)
