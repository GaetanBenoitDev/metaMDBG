project(metaMDBG)
#project(metaMDBG VERSION 0.1)

cmake_minimum_required(VERSION 3.10)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -w") #-g 
set(CMAKE_CXX_FLAGS "-Wall -w -g -Wno-c++11-narrowing -D_LIBCPP_DISABLE_AVAILABILITY")

#set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(ZLIB)

set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*)
add_executable(${PROJECT_NAME} ${ProjectFiles})
#target_link_libraries(${PROJECT_NAME})

#set (PROGRAM_SOURCE_DIR_MAPPER ${PROJECT_SOURCE_DIR}/src_mapper)
#include_directories (${PROGRAM_SOURCE_DIR_MAPPER})
#file (GLOB_RECURSE  ProjectFilesMapper  ${PROGRAM_SOURCE_DIR_MAPPER}/*)
#add_executable(mapper ${ProjectFilesMapper})



#if (NOT TARGET spoa)
#  add_subdirectory(src/utils/spoa64)
#endif ()
#target_link_libraries(${PROJECT_NAME} spoa64::spoa64)


if (NOT TARGET spoa)
  add_subdirectory(src/utils/spoa)
endif ()
target_link_libraries(${PROJECT_NAME} spoa::spoa)

#if (NOT TARGET abPOA)
#  add_subdirectory(src/utils/abPOA2 EXCLUDE_FROM_ALL)
#endif ()
#target_link_libraries(${PROJECT_NAME} abpoa)

target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
#target_link_libraries(mapper ZLIB::ZLIB)

