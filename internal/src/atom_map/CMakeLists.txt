cmake_minimum_required(VERSION 2.8.3)
project(atom_map)

find_package(PCL REQUIRED COMPONENTS common)
if(NOT PCL_FOUND)
  message(FATAL_ERROR "This program requires the PCL library.")
endif(NOT PCL_FOUND)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  kdtree
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
    roscpp
    kdtree
  DEPENDS
    PCL
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)

link_directories(
  ${catkin_LIBRARY_DIRS}
  ${PCL_LIBRARY_DIRS}
)

file(GLOB_RECURSE ${PROJECT_NAME}_srcs ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src/*.cc)
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_srcs})
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
)

if(CATKIN_ENABLE_TESTING)
  catkin_add_gtest(test_${PROJECT_NAME} test/test_atom_map.cc)
  target_link_libraries(test_${PROJECT_NAME} ${PROJECT_NAME} ${GTEST_LIBRARIES})
endif()
