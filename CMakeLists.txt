# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.2...3.12)
project (example CXX)

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

# ----------------------------------------------------------------------------
# Define paths
# ----------------------------------------------------------------------------

set(SEQAN_ROOT "${CMAKE_SOURCE_DIR}/submodules/seqan/")
set(SDSL_ROOT "${CMAKE_SOURCE_DIR}/submodules/sdsl-lite/")
set(BENCHMARK_ROOT "${CMAKE_SOURCE_DIR}/submodules/benchmark/")

# SeqAn include
set (CMAKE_INCLUDE_PATH ${SEQAN_ROOT}/include ${CMAKE_INCLUDE_PATH})

# SeqAn cmake module
set (CMAKE_PREFIX_PATH ${SEQAN_ROOT}/util/cmake ${CMAKE_PREFIX_PATH})
set (CMAKE_MODULE_PATH ${SEQAN_ROOT}/util/cmake ${CMAKE_MODULE_PATH})

include_directories (${CMAKE_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# Require C++17
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_ORIGINAL ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
    "#if !defined (__cplusplus) || (__cplusplus < 201703L)
    #error NOCXX17
    #endif
    int main() {}")

check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX17_BUILTIN)

if (CXX17_BUILTIN)
    message (STATUS "C++ Standard-17 support:    builtin")
else ()
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_ORIGINAL} -std=c++17")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX17_FLAG)

    if (CXX17_FLAG)
        message (STATUS "C++ Standard-17 support:    via -std=c++17")
    else ()
        message (FATAL_ERROR "SeqAn3 requires C++17, but your compiler does not support it.")
    endif ()

    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif ()

# ----------------------------------------------------------------------------
# Find zlib
# ----------------------------------------------------------------------------

find_package (ZLIB REQUIRED)

# ----------------------------------------------------------------------------
# Load SeqAn module
# ----------------------------------------------------------------------------

find_package (SeqAn REQUIRED)

# add include dir of test directory for some global test utility functions.
set (SEQAN_INCLUDE_DIRS "${SEQAN_INCLUDE_DIRS};${CMAKE_SOURCE_DIR}/include")

# import package's information (libraries are added with targets below)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

# add warning levels
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror -Wall -Wextra")

# ----------------------------------------------------------------------------
# Use CCache if possible
# ----------------------------------------------------------------------------

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif(CCACHE_PROGRAM)

# ----------------------------------------------------------------------------
# Add google/benchmark
# ----------------------------------------------------------------------------

set(google_benchmark_args "")
list (APPEND google_benchmark_args "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND google_benchmark_args "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND google_benchmark_args "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND google_benchmark_args "-DBENCHMARK_ENABLE_TESTING=false")
list (APPEND google_benchmark_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

include(ExternalProject)
ExternalProject_Add(
    google_benchmark
    PREFIX google_benchmark
    SOURCE_DIR "${BENCHMARK_ROOT}"
    CMAKE_ARGS "${google_benchmark_args}"
)

unset (google_benchmark_args)

add_library (benchmark STATIC IMPORTED)
add_dependencies(benchmark google_benchmark)
set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libbenchmark.a")

# ----------------------------------------------------------------------------
# define all include folders
# ----------------------------------------------------------------------------

include_directories("${BENCHMARK_ROOT}/include/")
include_directories("${SEQAN_ROOT}/include/")
include_directories("${SDSL_ROOT}/include/")
include_directories("${CMAKE_SOURCE_DIR}/include/")

# ----------------------------------------------------------------------------
# Add the executables
# ----------------------------------------------------------------------------

add_definitions(-DBASE_DIR="${CMAKE_SOURCE_DIR}")

add_executable(minimizer_example src/minimizer_example.cpp)
target_link_libraries(minimizer_example "pthread")

add_executable(minimizer_benchmark src/minimizer_benchmark.cpp)
target_link_libraries(minimizer_benchmark "pthread" "benchmark" "z")
add_dependencies(minimizer_benchmark benchmark)

add_executable(ibf_example src/ibf_example.cpp)
target_link_libraries(ibf_example "pthread")
