cmake_minimum_required(VERSION 2.6.2)
project(tgcq CXX C)

### Flags for all
SET(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

# Where .h files can be found.
#include_directories(
#  ${CMAKE_CURRENT_SOURCE_DIR}/include
#  ${CMAKE_CURRENT_SOURCE_DIR}/src)

# libcqtree
add_subdirectory(../libcqtree libcqtree_build)
include_directories(${libcqtree_SOURCE_DIR}/src)

# libcds used to build libcqtree
include_directories(${libcqtree_SOURCE_DIR}/libcds_ext-prefix/include)
add_library(libcds STATIC IMPORTED)
SET_TARGET_PROPERTIES(libcds PROPERTIES
    IMPORTED_LOCATION  ${CMAKE_CURRENT_SOURCE_DIR}/libcqtree_build/libcds_ext-prefix/lib/libcds.a)

add_dependencies(libcds libcds_ext)



add_executable( create create.cc $<TARGET_OBJECTS:Objs>)
target_link_libraries(create libcds)

### Common objects files in my project
#add_library(Objs OBJECT src/utils.cc src/MXCompactQtree.cc src/MXCompactQtreeFixed.cc src/PRBCompactQtree.cc src/PRWCompactQtree.cc)