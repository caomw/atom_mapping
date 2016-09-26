cmake_minimum_required(VERSION 2.8.3)
project(atom_map_example)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(PCL REQUIRED COMPONENTS common filters)
if(NOT PCL_FOUND)
  message(FATAL_ERROR "This program requires the PCL library.")
endif(NOT PCL_FOUND)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  pcl_ros
  tf2_ros
  atom_map
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
    roscpp
    pcl_ros
    tf2_ros
    atom_map
  DEPENDS
    PCL
    FLANN
)

option(ENABLE_DEBUG_MESSAGES "Turn on to enable bound checking and debug messages" OFF)
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

add_library(${PROJECT_NAME} src/${PROJECT_NAME}.cc)
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${FLANN_LIBRARIES}
  glog
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
)

add_executable(${PROJECT_NAME}_node src/${PROJECT_NAME}_node.cc)
target_link_libraries(${PROJECT_NAME}_node
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
)
