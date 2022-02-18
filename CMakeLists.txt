project(mdbgAsmMeta)
#project(mdbgAsmMeta VERSION 0.1)

cmake_minimum_required(VERSION 3.10)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(ZLIB)

set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*)
add_executable(${PROJECT_NAME} ${ProjectFiles})
#target_link_libraries(${PROJECT_NAME})


if (NOT TARGET spoa)
  add_subdirectory(src/utils/spoa EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(${PROJECT_NAME} spoa::spoa)
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)

