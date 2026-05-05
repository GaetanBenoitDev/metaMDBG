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

set(TURBOPFOR_DIR ${PROJECT_SOURCE_DIR}/ext/TurboPFor-Integer-Compression)
set(TURBOPFOR_LIB ${TURBOPFOR_DIR}/libic.a)

execute_process(
    COMMAND make -j 4
    WORKING_DIRECTORY ${TURBOPFOR_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "TurboPFor build failed")
endif()



add_custom_target(turbopfor_build
    DEPENDS ${TURBOPFOR_LIB}
)

add_library(turbopfor_lib STATIC IMPORTED)

set_target_properties(turbopfor_lib PROPERTIES
    IMPORTED_LOCATION ${TURBOPFOR_LIB}
)

target_link_libraries(${PROJECT_NAME} turbopfor_lib)





















set(MINIMAP2_DIR ${PROJECT_SOURCE_DIR}/ext/minimap2)
set(MINIMAP2_LIB ${MINIMAP2_DIR}/libminimap2.a)

execute_process(
    COMMAND make -j 4
    WORKING_DIRECTORY ${MINIMAP2_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "Minimap2 build failed")
endif()

add_custom_target(minimap2_build
    DEPENDS ${MINIMAP2_LIB}
)

add_library(minimap2_lib STATIC IMPORTED)

set_target_properties(minimap2_lib PROPERTIES
    IMPORTED_LOCATION ${MINIMAP2_LIB}
)

target_link_libraries(${PROJECT_NAME} minimap2_lib)












set(SPOA_DIR ${PROJECT_SOURCE_DIR}/ext/spoa)
set(SPOA_LIB ${SPOA_DIR}/build/lib/libspoa.a)

execute_process(
    COMMAND cmake -B build -DCMAKE_BUILD_TYPE=Release
    WORKING_DIRECTORY ${SPOA_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "Spoa build failed")
endif()

execute_process(
    COMMAND make -C build -j 4
    WORKING_DIRECTORY ${SPOA_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "Spoa build failed")
endif()

add_custom_target(spoa_build
    DEPENDS ${SPOA_LIB}
)

add_library(spoa_lib STATIC IMPORTED)

set_target_properties(spoa_lib PROPERTIES
    IMPORTED_LOCATION ${SPOA_LIB}
)

target_link_libraries(${PROJECT_NAME} spoa_lib)






set(HTSLIB_DIR ${PROJECT_SOURCE_DIR}/ext/htslib)
set(HTSLIB_LIB ${HTSLIB_DIR}/libhts.a)

execute_process(
    COMMAND make -j 4
    WORKING_DIRECTORY ${HTSLIB_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "Htslib build failed")
endif()

add_custom_target(htslib_build
    DEPENDS ${HTSLIB_LIB}
)

add_library(htslib_lib STATIC IMPORTED)

set_target_properties(htslib_lib PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_LIB}
)

target_link_libraries(${PROJECT_NAME} htslib_lib)





#if (NOT TARGET spoa)
#  add_subdirectory(src/utils/spoa)
#endif ()
#target_link_libraries(${PROJECT_NAME} spoa::spoa)




#if (NOT TARGET abPOA)
#  add_subdirectory(src/utils/abPOA EXCLUDE_FROM_ALL)
#endif ()
#target_link_libraries(${PROJECT_NAME} abpoa)

target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
#target_link_libraries(mapper ZLIB::ZLIB)





