cmake_minimum_required(VERSION 2.6.2)
project(tgcq CXX C)

### Flags for all
SET(CMAKE_CXX_FLAGS "-std=c++11  -O3 -Wall -Wextra")

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
# this include sdsl cds cdswrapper
add_subdirectory(libcqtree)

###
# cdswrapper
##
#add_subdirectory(cdswrapper)
#include_directories(cdswrapper/include)

##
#libcds
##
#add_subdirectory(../libcds ext_libcds)
#include_directories(../libcds/include)

##
#sdsl
##
#add_subdirectory(../sdsl-lite ext_sdsl-lite)
#include_directories(../sdsl-lite/include)

add_executable( create create.cc )
target_link_libraries(create cqtree)

add_executable( sort sort.cc )
target_link_libraries(sort cqtree)


add_executable( benchmark benchmark.cc arraysort.cc TemporalGraph.cc)
target_link_libraries(benchmark cqtree)
set_target_properties (benchmark PROPERTIES COMPILE_DEFINITIONS "EXPERIMENTS")

add_executable( checkans benchmark.cc arraysort.cc TemporalGraph.cc)
target_link_libraries(checkans cqtree)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        target_link_libraries(benchmark -lrt)
	target_link_libraries(checkans -lrt)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
