cmake_minimum_required(VERSION 2.8.3)
project(atom_map)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PCL REQUIRED COMPONENTS common features)
if(NOT PCL_FOUND)
  message(FATAL_ERROR "This program requires the PCL library.")
endif(NOT PCL_FOUND)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  pcl_ros
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
    roscpp
    pcl_ros
  DEPENDS
    PCL
    FLANN
)

option(ENABLE_DEBUG_MESSAGES "Turn on to enable bound checking and debug messages" ON)
if(ENABLE_DEBUG_MESSAGES)
  add_definitions(-DENABLE_DEBUG_MESSAGES=1)
else()
  add_definitions(-DENABLE_DEBUG_MESSAGES=0)
endif()

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
  ${FLANN_INCLUDE_DIRS}
  ${GLOG_INCLUDE_DIRS}
)

link_directories(
  ${catkin_LIBRARY_DIRS}
  ${PCL_LIBRARY_DIRS}
  ${FLANN_LIBRARY_DIRS}
  ${GLOG_LIBRARY_DIRS}
)

file(GLOB_RECURSE ${PROJECT_NAME}_srcs ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src/*.cc)
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_srcs})
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${FLANN_LIBRARIES}
  glog
)

if(CATKIN_ENABLE_TESTING)
  catkin_add_gtest(test_${PROJECT_NAME} test/test_atom_map.cc)
  target_link_libraries(
    test_${PROJECT_NAME}
    ${PROJECT_NAME}
    ${PCL_LIBRARIES}
    ${FLANN_LIBRARIES}
    ${GTEST_LIBRARIES}
  )
endif()
