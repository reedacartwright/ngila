###########################################################
#  Ngila's Project Configuration
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

set(NEW_PACKAGE_NAME "ngila")
set(NEW_PACKAGE_VERSION "1-current-rUnknown")

SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Turn on Verbose Makefiles" FORCE)
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

PROJECT(Ngila)

#### Static Lib Support ####

SET(USE_STATIC_LIBS OFF CACHE BOOL "Link against static libraries.")
IF(USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_USE_STATIC_RUNTIME ON)
ENDIF(USE_STATIC_LIBS)

#SET(USE_THREADS OFF CACHE BOOL
#	"Enable multi-threading support in Ngila (experimental).")
#MARK_AS_ADVANCED(USE_THREADS)

#### RPATH Support #####

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()


INCLUDE(NewPackage)

if(WIN32 AND NOT UNIX)
	set(DIR_SUFFIX "")
else()
	set(DIR_SUFFIX "\${NEW_PACKAGE_INSTALL_SUFFIX}")
endif()
SET(CMAKE_DATA_DIR "share${DIR_SUFFIX}" CACHE STRING
  "Install location for data (relative to prefix).")
SET(CMAKE_DOC_DIR "doc${DIR_SUFFIX}" CACHE STRING
  "Install location for documentation (relative to prefix).")
SET(CMAKE_MAN_DIR "man" CACHE STRING
  "Install location for man pages (relative to prefix).")
MARK_AS_ADVANCED(CMAKE_DATA_DIR CMAKE_DOC_DIR CMAKE_MAN_DIR)

SET(CPACK_SOURCE_IGNORE_FILES
  "ngila\\\\.exe"  "ngila$"
)

SET(CPACK_PACKAGE_NAME "Ngila")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pairwise alignment with evolutionary models.")
SET(CPACK_PACKAGE_VENDOR "Reed A. Cartwright")
SET(CPACK_PACKAGE_EXECUTABLES "ngila" "Ngila")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})

IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  SET(CPACK_NSIS_DISPLAY_NAME "Ngila ${NEW_PACKAGE_VERSION}")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\ngila.exe")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\scit.us\\\\projects\\\\ngila\\\\")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\scit.us\\\\projects\\\\ngila\\\\")
  SET(CPACK_NSIS_CONTACT "cartwright@asu.edu")
ELSE()
  SET(CPACK_STRIP_FILES "bin/ngila")
ENDIF()

INCLUDE(CPackBefore)
INCLUDE(CPack)

INSTALL(FILES copying.txt readme.txt ngilarc.txt ChangeLog DESTINATION ${CMAKE_DOC_DIR})
#INSTALL(FILES lambda.pl outsplit.pl nexus2fasta.pl varrep.pl DESTINATION ${CMAKE_DATA_DIR})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(matrix)
ADD_SUBDIRECTORY(doc)
#ADD_SUBDIRECTORY(tests)

ADD_DEPENDENCIES(new_package ngila)
ADD_DEPENDENCIES(ngila version)

########### Add uninstall target ###############
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/Modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake") 

#INCLUDE(InstallRequiredSystemLibraries)


