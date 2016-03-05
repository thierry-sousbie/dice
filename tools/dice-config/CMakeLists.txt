set(MAIN_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/")
set(MAIN_BINARY_NAME "dice-config")

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY_BACK ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

SET( DICECONFIG_SRC
  dice-config.cpp  
  dice-config.hxx 
  )

SET( MAIN_BINARY_SRC
  ${DICECONFIG_SRC}
  )

add_executable("${MAIN_BINARY_NAME}" ${MAIN_BINARY_SRC})

file(WRITE dice-config.hxx "char CXX_FLAGS[]=\"${CMAKE_CXX_FLAGS}\";\n" )
#file(APPEND dice-config.hxx "char LINKER_FLAGS[]=\"${CMAKE_EXE_LINKER_FLAGS}\";\n" )

get_property(myDefinesVar 
  DIRECTORY 
  PROPERTY COMPILE_DEFINITIONS)
string(REPLACE ";" " -D" myDefines "${myDefinesVar}")
file(APPEND dice-config.hxx "char DEFINES[]=\" -D${myDefines}\";\n" )

get_property(myIncludesVar 
  DIRECTORY 
  PROPERTY INCLUDE_DIRECTORIES)
string(REPLACE ";" " -I" myIncludes "${myIncludesVar}")
string(REPLACE "${PROJECT_SOURCE_DIR}" "${CMAKE_INSTALL_PREFIX}" myInstalledIncludes "${myIncludes}")
string(REPLACE "/external" "/dice/external" myInstalledIncludes "${myInstalledIncludes}")

file(APPEND dice-config.hxx "char INCLUDES[]=\" -I${myInstalledIncludes}\";\n" )
file(APPEND dice-config.hxx "char DIST_INCLUDES[]=\" -I${myIncludes}\";\n" )

get_property(myLibrariesVar 
  TARGET ${MAIN_BINARY_NAME}
  PROPERTY LINK_LIBRARIES)
string(REPLACE ";" " " myLibraries "${myLibrariesVar}")

file(APPEND dice-config.hxx "char LIBRARIES[]=\" ${myLibraries}\";\n" )

file(APPEND dice-config.hxx "char COMPILER[]=\" ${CMAKE_CXX_COMPILER}\";\n" )

set(ALL_TARGETS "${MAIN_BINARY_NAME}")

SET(TARGETS_LIST_LOC ${TARGETS_LIST_LOC} ${ALL_TARGETS})

install (TARGETS ${ALL_TARGETS} DESTINATION bin)

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(TARGETS_LIST ${TARGETS_LIST} ${TARGETS_LIST_LOC} PARENT_SCOPE)
SET(NOT_BUILT_TARGETS_LIST ${NOT_BUILT_TARGETS_LIST} ${NOT_BUILT_TARGETS_LIST_LOC} PARENT_SCOPE)
SET(DISABLED_FEATURES_LIST ${DISABLED_FEATURES_LIST} ${DISABLED_FEATURES_LIST_LOC} PARENT_SCOPE)
