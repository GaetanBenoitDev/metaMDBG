cmake_minimum_required(VERSION 3.2)
project(abpoa LANGUAGES C VERSION 3.0.0)

include(GNUInstallDirs)
find_package(ZLIB REQUIRED)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native") # SIMD

# build abPOA as a static library by default
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

add_library(abpoa
    src/abpoa_align.c
    src/abpoa_graph.c
    src/abpoa_output.c
    src/abpoa_plot.c
    src/abpoa_seed.c
    src/abpoa_seq.c
    src/kalloc.c
    src/kstring.c
    src/simd_abpoa_align.c
    src/simd_check.c
    src/utils.c)


add_executable(abpoa_bin
    src/abpoa.c)

target_link_libraries(abpoa_bin abpoa z pthread m)
set_target_properties(abpoa_bin PROPERTIES OUTPUT_NAME abpoa)

target_include_directories(abpoa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:src>)


install(TARGETS abpoa DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS abpoa_bin DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")

# configure and install pkg-config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/abpoa.pc.in ${CMAKE_CURRENT_BINARY_DIR}/abpoa-1.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/abpoa-1.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
