cmake_minimum_required(VERSION 3.1)
project(reTERRA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER:FILEPATH C:/msys64/mingw64/bin/g++.exe)
set(CMAKE_C_COMPILER:FILEPATH C:/msys64/mingw64/bin/gcc.exe)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


#This project requires CGAL
find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()
include_directories(${CGAL_INCLUDE_DIRS})


#CGAL requires Boost
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

#This project requires Eigen
find_package( Eigen3 3.4.0 REQUIRED NO_MODULE)
if ( NOT EIGEN3_FOUND )
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()
include_directories(${Eigen3_INCLUDE_DIRS})
message(STATUS "Found Eigen3 include dir:" ${EIGEN3_INCLUDE_DIRS})

find_package(OpenCV REQUIRED)
if ( NOT OpenCV_FOUND )
  message(STATUS "This project requires the OpenCV library, and will not be compiled.")
  return()
endif()
include_directories(${OpenCV_INCLUDE_DIRS})
message(STATUS "Found OpenCV include dir:" ${OpenCV_INCLUDE_DIRS})
message(STATUS "Found OpenCV lib:" ${OpenCV_LIBS})

# Add this line if liblkhUtils.a is missing or broken
# add_subdirectory(lkhUtils lkhUtils/OBJ)

link_directories(${PROJECT_SOURCE_DIR}/lkhUtils/OBJ)

include_directories(
    ${PROJECT_SOURCE_DIR}/include

)
add_executable(${PROJECT_NAME}
    ${PROJECT_SOURCE_DIR}/src/TERRAExec.cpp
    ${PROJECT_SOURCE_DIR}/src/TERRA.cpp
    ${PROJECT_SOURCE_DIR}/src/Voronoi.cpp
    ${PROJECT_SOURCE_DIR}/src/SetCover.cpp
    ${PROJECT_SOURCE_DIR}/src/TspGaUgv.cpp
    ${PROJECT_SOURCE_DIR}/src/UAVComputePath.cpp
    ${PROJECT_SOURCE_DIR}/src/SearchUAVOperations.cpp
    ${PROJECT_SOURCE_DIR}/src/TSP.cpp
    ${PROJECT_SOURCE_DIR}/src/BuildSolution.cpp
)

target_link_libraries(${PROJECT_NAME} lkhUtils)
target_link_libraries(${PROJECT_NAME} ${CGAL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})