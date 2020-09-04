# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.2...3.12)
project (benchmark CXX)

# ----------------------------------------------------------------------------
# Find zlib
# ----------------------------------------------------------------------------

find_package (ZLIB REQUIRED)

# ----------------------------------------------------------------------------
# Find SeqAn3
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED HINTS ${CMAKE_SOURCE_DIR}/submodules/seqan3/build_system/)

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

# set(google_benchmark_args "")
# list (APPEND google_benchmark_args "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
# list (APPEND google_benchmark_args "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
# list (APPEND google_benchmark_args "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
# list (APPEND google_benchmark_args "-DBENCHMARK_ENABLE_TESTING=false")
# list (APPEND google_benchmark_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")
#
# include(ExternalProject)
# ExternalProject_Add(
#     google_benchmark
#     PREFIX google_benchmark
#     SOURCE_DIR "${BENCHMARK_ROOT}"
#     CMAKE_ARGS "${google_benchmark_args}"
# )
#
# unset (google_benchmark_args)
#
# add_library (benchmark STATIC IMPORTED)
# add_dependencies(benchmark google_benchmark)
# set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libbenchmark.a")

# ----------------------------------------------------------------------------
# define all include folders
# ----------------------------------------------------------------------------

include_directories("${CMAKE_SOURCE_DIR}/include/")
# include_directories("${BENCHMARK_ROOT}/include/")

# ----------------------------------------------------------------------------
# Add the executables
# ----------------------------------------------------------------------------

# add_definitions(-DBASE_DIR="${CMAKE_SOURCE_DIR}")

# add_executable(minimizer_example src/minimizer_example.cpp)
# target_link_libraries(minimizer_example "pthread" seqan3::seqan3)

# add_executable(minimizer_benchmark src/minimizer_benchmark.cpp)
# target_link_libraries(minimizer_benchmark "pthread" "benchmark" "z")
# add_dependencies(minimizer_benchmark benchmark seqan3::seqan3)

# add_executable(ibf_example src/ibf_example.cpp)
# target_link_libraries(ibf_example "pthread" seqan3::seqan3)

add_executable(build_min_ibf src/build_min_ibf.cpp)
target_link_libraries(build_min_ibf seqan3::seqan3)

add_executable(search_min_ibf src/search_min_ibf.cpp)
target_link_libraries(search_min_ibf seqan3::seqan3)

add_executable(generate_reads src/generate_reads.cpp)
target_link_libraries(generate_reads seqan3::seqan3)
