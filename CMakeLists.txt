# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(odgi)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -DGAINBUCKET_ARRAY -DMONDRIAANVERSION='\"4.2.1\"'")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DGAINBUCKET_ARRAY -DMONDRIAANVERSION='\"4.2.1\"'")

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_library(mondriaan_objs OBJECT
  ${CMAKE_SOURCE_DIR}/src/Cartesian.c
  ${CMAKE_SOURCE_DIR}/src/DistributeMat.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVec.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecGreedy.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecLib.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecLocal.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecOpt2.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecOrig.c
  ${CMAKE_SOURCE_DIR}/src/DistributeVecOrigEq.c
  ${CMAKE_SOURCE_DIR}/src/FreeNonzeros.c
  ${CMAKE_SOURCE_DIR}/src/GainBucketArray.c
  ${CMAKE_SOURCE_DIR}/src/Graph.c
  ${CMAKE_SOURCE_DIR}/src/Heap.c
  ${CMAKE_SOURCE_DIR}/src/HKLFM.c
  ${CMAKE_SOURCE_DIR}/src/Matalloc.c
  ${CMAKE_SOURCE_DIR}/src/Match.c
  ${CMAKE_SOURCE_DIR}/src/MatchInproduct.c
  ${CMAKE_SOURCE_DIR}/src/MatchMatchers.c
  ${CMAKE_SOURCE_DIR}/src/MatchStairway.c
  ${CMAKE_SOURCE_DIR}/src/Options.c
  ${CMAKE_SOURCE_DIR}/src/Permute.c
  ${CMAKE_SOURCE_DIR}/src/Remembrance.c
  ${CMAKE_SOURCE_DIR}/src/Sort.c
  ${CMAKE_SOURCE_DIR}/src/SparseMatrix.c
  ${CMAKE_SOURCE_DIR}/src/SplitMatrixUpperBound.c
  ${CMAKE_SOURCE_DIR}/src/SubsetSum.c
  ${CMAKE_SOURCE_DIR}/src/ZeroVolumeSearch.c)
add_library(libmondriaan STATIC $<TARGET_OBJECTS:mondriaan_objs>)
set_target_properties(libmondriaan PROPERTIES OUTPUT_NAME "mondriaan")

add_executable(mondriaan
  $<TARGET_OBJECTS:mondriaan_objs>
  ${CMAKE_SOURCE_DIR}/tools/Mondriaan.c
  )
target_include_directories(mondriaan PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}/tools")
target_link_libraries(mondriaan
  -lm
  )

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
