cmake_minimum_required(VERSION 2.6.2)
project(tgcq CXX C)

### Flags for all
SET(CMAKE_CXX_FLAGS "")

# Where .h files can be found.
#include_directories(
#  ${CMAKE_CURRENT_SOURCE_DIR}/include
#  ${CMAKE_CURRENT_SOURCE_DIR}/src)

### Default parameter for debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O2 -Wall -Wextra")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   message(STATUS "Makefile for Debug mode")
endif()

### Default parametrs for release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
   message(STATUS "Makefile for Release mode")
endif()



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

add_executable( sort sort.cc )

add_definitions( -DEXPERIMENTS )

add_executable( benchmark benchmark.cc arraysort.cc $<TARGET_OBJECTS:Objs>)
target_link_libraries(benchmark libcds)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(benchmark -lrt)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

### Common objects files in my project
#add_library(Objs OBJECT src/utils.cc src/MXCompactQtree.cc src/MXCompactQtreeFixed.cc src/PRBCompactQtree.cc src/PRWCompactQtree.cc)
