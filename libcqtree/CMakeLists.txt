cmake_minimum_required(VERSION 2.6.2)
project(libcqtree CXX C)

### Flags for all
SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

# Where .h files can be found.
include_directories(
  ${libcqtree_SOURCE_DIR}/include
  ${libcqtree_SOURCE_DIR}/src)

# Where libraries can be found.
# link_directories(${gtest_BINARY_DIR}/src)

################################
# GTest
################################
add_subdirectory (../external/gtest-1.7.0 ext_gtest-1.7.0)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

###
# cdswrapper, sdsl and libcds are automatically included :D
##
add_subdirectory(../cdswrapper ext_cdswrapper)
#include_directories(../cdswrapper/include)

#libcds
#add_subdirectory(../../libcds ext_libcds)
#include_directories(../../libcds/include)

#sdsl
#add_subdirectory(../../sdsl-lite ext_sdsl-lite)
#include_directories(../../sdsl-lite/include)

################################
# Unit Tests
################################
enable_testing()
# Add test cpp file

add_executable( runPointUnitTests test/point_test.cc )
# Link test executable against gtest & gtest_main
target_link_libraries(runPointUnitTests gtest gtest_main)
add_test( runPointUnitTests runPointUnitTests )

add_executable( runCodeUnitTests test/code_test.cc )
# Link test executable against gtest & gtest_main
target_link_libraries(runCodeUnitTests gtest gtest_main cqtree)
add_test( runCodeUnitTests runCodeUnitTests )

add_executable( runMXUnitTests test/mx_test.cc  )
target_link_libraries(runMXUnitTests gtest gtest_main cqtree)
add_test( runMXUnitTests runMXUnitTests )

#add_executable( runMXFUnitTests test/mxf_test.cc  )
#target_link_libraries(runMXFUnitTests gtest gtest_main cqtree)
#add_test( runMXFUnitTests runMXFUnitTests )

add_executable( runPRBUnitTests test/prb_test.cc  )
target_link_libraries(runPRBUnitTests gtest gtest_main cqtree)
add_test( runPRBUnitTests runPRBUnitTests )

#add_executable( runPRWUnitTests test/prw_test.cc  )
#target_link_libraries(runPRWUnitTests gtest gtest_main cqtree)
#add_test( runPRWUnitTests runPRWUnitTests )

### Common objects files in my project
#add_library(cqtree src/utils.cc src/CompactQtree.cc src/MXCompactQtree.cc src/MXCompactQtreeFixed.cc src/PRBCompactQtree.cc src/PRWCompactQtree.cc )
add_library(cqtree src/utils.cc src/CompactQtree.cc src/MXCompactQtree.cc src/PRBCompactQtree.cc)

target_include_directories(cqtree PUBLIC src)
target_include_directories(cqtree PUBLIC include)

target_link_libraries(cqtree cdswrapper)