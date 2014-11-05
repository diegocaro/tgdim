cmake_minimum_required(VERSION 2.6.2)
project(tgcq CXX C)

### Flags for all
SET(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++  -O3 -Wall -Wextra")

# Where .h files can be found.
#include_directories(
#  ${CMAKE_CURRENT_SOURCE_DIR}/include
#  ${CMAKE_CURRENT_SOURCE_DIR}/src)

### Default parameter for debug
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2 -Wall -Wextra")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   message(STATUS "Makefile for Debug mode")
   #    add_definitions( -DEXPERIMENTS )
endif()

### Default parametrs for release
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ")
#if (CMAKE_BUILD_TYPE STREQUAL "Release")
#   message(STATUS "Makefile for Release mode")
#endif()



# libcqtree
add_subdirectory(../libcqtree libcqtree_build)
include_directories(${libcqtree_SOURCE_DIR}/src)

# libcds used to build libcqtree
include_directories(${libcqtree_SOURCE_DIR}/libcds_ext-prefix/include)
add_library(libcds STATIC IMPORTED)
SET_TARGET_PROPERTIES(libcds PROPERTIES
    IMPORTED_LOCATION  ${CMAKE_CURRENT_SOURCE_DIR}/libcqtree_build/libcds_ext-prefix/lib/libcds.a)

add_dependencies(libcds libcds_ext)

# libsdsl used to build libcqtree
include_directories(${libcqtree_SOURCE_DIR}/../sdsl-lite/include)
add_library(libsdsl STATIC IMPORTED)
SET_TARGET_PROPERTIES(libsdsl PROPERTIES
    IMPORTED_LOCATION  ${CMAKE_CURRENT_SOURCE_DIR}/../sdsl-lite/lib/libsdsl.a)

add_dependencies(libsdsl libsdsl_ext)



add_executable( create create.cc $<TARGET_OBJECTS:Objs>)
target_link_libraries(create libcds libsdsl)

add_executable( sort sort.cc )



add_executable( benchmark benchmark.cc arraysort.cc TemporalGraph.cc $<TARGET_OBJECTS:Objs>)
target_link_libraries(benchmark libcds libsdsl)
set_target_properties (benchmark PROPERTIES COMPILE_DEFINITIONS "EXPERIMENTS")


add_executable( checkans benchmark.cc arraysort.cc TemporalGraph.cc $<TARGET_OBJECTS:Objs>)
target_link_libraries(checkans libcds libsdsl)


IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(benchmark -lrt)
		target_link_libraries(checkans -lrt)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

### Common objects files in my project
#add_library(Objs OBJECT src/utils.cc src/MXCompactQtree.cc src/MXCompactQtreeFixed.cc src/PRBCompactQtree.cc src/PRWCompactQtree.cc)
