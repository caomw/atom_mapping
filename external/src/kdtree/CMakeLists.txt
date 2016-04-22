cmake_minimum_required(VERSION 2.8)
project(kdtree)
find_package(catkin REQUIRED)
catkin_package()

include(ExternalProject)

set(kdtree_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/kdtree")
set(kdtree_SOURCE_DIR "${kdtree_PREFIX}/src/kdtree")
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${kdtree_PREFIX}/src)

externalproject_add(kdtree_download
  PREFIX ${kdtree_PREFIX}
  SOURCE_DIR ${kdtree_SOURCE_DIR}
  GIT_REPOSITORY https://github.com/jtsiomb/kdtree.git
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
)

set_source_files_properties(${kdtree_SOURCE_DIR}/kdtree.c PROPERTIES GENERATED 1)
set_source_files_properties(${kdtree_SOURCE_DIR}/kdtree.h PROPERTIES GENERATED 1)

add_library(kdtree SHARED ${kdtree_SOURCE_DIR}/kdtree.c)
add_dependencies(kdtree kdtree_download)

set(kdtree_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
set(kdtree_LIBRARY_DIRS "${CMAKE_INSTALL_PREFIX}/lib")

install(FILES ${kdtree_SOURCE_DIR}/kdtree.h DESTINATION "${kdtree_INCLUDE_DIRS}/kdtree")
install(TARGETS kdtree DESTINATION ${kdtree_LIBRARY_DIRS})

install(FILES package.xml DESTINATION share/${PROJECT_NAME})
