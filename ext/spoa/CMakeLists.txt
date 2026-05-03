cmake_minimum_required(VERSION 3.12)

project(spoa VERSION 4.1.5
             LANGUAGES CXX
             DESCRIPTION "Spoa is a c++ library (and tool) for SIMD vectorized partial order alignment.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(spoa_main_project ON)
endif ()
option(spoa_install "Generate library install target" ${spoa_main_project})
option(spoa_build_exe "Build executable" ${spoa_main_project})
option(spoa_build_tests "Build unit tests" ${spoa_main_project})
option(spoa_use_cereal "Use cereal library" OFF)
option(spoa_optimize_for_native "Build with -march=native" ON)
option(spoa_optimize_for_portability "Build with -msse4.1" OFF)
option(spoa_use_simde "Use SIMDe library for porting vectorized code" OFF)
option(spoa_use_simde_nonvec "Use SIMDe library for nonvectorized code" OFF)
option(spoa_use_simde_openmp "Use SIMDe support for OpenMP SIMD" OFF)
option(spoa_generate_dispatch "Use SIMDe to generate x86 dispatch" OFF)

if (NOT spoa_generate_dispatch)
  if (spoa_optimize_for_portability)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
  elseif (spoa_optimize_for_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif ()
endif ()
if (spoa_use_simde OR
    spoa_use_simde_nonvec OR
    spoa_use_simde_openmp OR
    spoa_generate_dispatch)
  FetchContent_Declare(
    simde
    GIT_REPOSITORY https://github.com/simd-everywhere/simde
    GIT_TAG v0.7.6)

  FetchContent_GetProperties(simde)
  if (NOT simde_POPULATED)
    FetchContent_Populate(simde)
  endif ()

  add_compile_definitions(SPOA_USE_SIMDE SIMDE_ENABLE_NATIVE_ALIASES)

  if (spoa_use_simde_nonvec)
    add_compile_definitions(SIMDE_NO_NATIVE)
  endif ()
  if (spoa_use_simde_openmp)
    add_compile_definitions(SIMDE_ENABLE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd")
  endif ()
  if (spoa_generate_dispatch)
    find_package(CpuFeatures 0.6.0 QUIET)

    if (NOT CpuFeatures_FOUND)
      if (spoa_install)
        message(FATAL_ERROR "Missing package CpuFeatures!")
      endif ()

      FetchContent_Declare(
        cpu_features
        GIT_REPOSITORY https://github.com/google/cpu_features
        GIT_TAG v0.6.0)

      FetchContent_GetProperties(cpu_features)
      if (NOT cpu_features_POPULATED)
        FetchContent_Populate(cpu_features)
        add_subdirectory(
          ${cpu_features_SOURCE_DIR}
          ${cpu_features_BINARY_DIR}
          EXCLUDE_FROM_ALL)
        add_library(CpuFeatures::cpu_features ALIAS cpu_features)
      endif ()
    endif ()

    add_compile_definitions(SPOA_GENERATE_DISPATCH)
  endif ()
endif ()

if (spoa_use_cereal)
  find_package(cereal 1.3.0 QUIET)

  if (NOT cereal_FOUND)
    if (spoa_install)
      message(FATAL_ERROR "Missing package cereal!")
    endif ()

    FetchContent_Declare(
      cereal
      GIT_REPOSITORY https://github.com/USCiLab/cereal
      GIT_TAG v1.3.0)

    FetchContent_GetProperties(cereal)
    if (NOT cereal_POPULATED)
      FetchContent_Populate(cereal)
      add_subdirectory(
        ${cereal_SOURCE_DIR}
        ${cereal_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(cereal::cereal ALIAS cereal)
    endif ()
  endif ()
endif ()


set(SPOA_VERSION "${PROJECT_VERSION}")
configure_file(src/spoa_config.h.in spoa_config.h)

add_library(spoa
  src/alignment_engine.cpp
  src/graph.cpp
  src/simd_alignment_engine_dispatcher.cpp
  src/sisd_alignment_engine.cpp
  src/version.cpp)
add_library(spoa::spoa ALIAS spoa)

target_include_directories(spoa PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${simde_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

if (spoa_use_cereal)
  target_link_libraries(spoa
    cereal::cereal)

  target_compile_definitions(spoa PUBLIC SPOA_USE_CEREAL)
endif ()

if (BUILD_SHARED_LIBS)
  set_property(TARGET spoa PROPERTY SOVERSION "7.0.0")
endif ()

if (spoa_generate_dispatch)
  list(APPEND ARCHITECTURES avx2 sse4.1 sse2)
  foreach(arch IN LISTS ARCHITECTURES)
    add_library(spoa_${arch} OBJECT
      src/simd_alignment_engine_dispatch.cpp)

    target_include_directories(spoa_${arch} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${simde_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>)

    set_property(TARGET spoa_${arch}
      PROPERTY COMPILE_FLAGS "-m${arch}")
    if (BUILD_SHARED_LIBS)
      set_property(TARGET spoa_${arch}
        PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif ()
  endforeach ()

  add_dependencies(spoa
    spoa_avx2
    spoa_sse4.1
    spoa_sse2)

  target_link_libraries(spoa
    CpuFeatures::cpu_features)
endif ()

if (spoa_install)
  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS spoa
    EXPORT ${PROJECT_NAME}Targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(
    DIRECTORY include/spoa
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif ()


