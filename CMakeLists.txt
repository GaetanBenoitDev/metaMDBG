project(Bloocoo)

cmake_minimum_required (VERSION 2.6)

################################################################################
# The version number.
################################################################################
# The default version number is the latest official build
SET (gatb-tool_VERSION_MAJOR 1)
SET (gatb-tool_VERSION_MINOR 0)
SET (gatb-tool_VERSION_PATCH 7)

# But, it is possible to define another release number during a local build
IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOLS
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories (${PROGRAM_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*)
add_executable(${PROJECT_NAME} ${ProjectFiles})
target_link_libraries(${PROJECT_NAME} ${gatb-core-libraries})
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ln -sf "${CMAKE_CURRENT_SOURCE_DIR}/scripts/metaBloocoo" "${PROJECT_BINARY_DIR}/bin"
)

################################################################################
#  PACKAGING
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# We ignore unwated files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES
    "^${PROJECT_SOURCE_DIR}/\\.git/"     ;
    "^${PROJECT_SOURCE_DIR}/\\.gitmodules"  ;
    "^${PROJECT_SOURCE_DIR}/\\.gitignore";
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${GATB_CORE_HOME}/\\.cproject" ;
    "^${GATB_CORE_HOME}/\\.git/"     ;
    "^${GATB_CORE_HOME}/\\.project"  ;
    "^${GATB_CORE_HOME}/\\.gitignore";
    "^${PROJECT_SOURCE_DIR}/tests/"  ;
    "^${PROJECT_SOURCE_DIR}/sam_tests/"
)

# We copy some files to the binary archive
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
INSTALL (DIRECTORY "${PROJECT_SOURCE_DIR}/test/" DESTINATION test)
INSTALL (FILES scripts/metaBloocoo DESTINATION bin)
INSTALL (FILES INSTALL LICENSE README.md DESTINATION bin/..)

# We include the "bin" tag into binary archive file name
set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME})

include (CPack)
