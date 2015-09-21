# Try to find SparseHash
# Once done, this will define
#
# SPARSEHASH_FOUND - system has SparseHash
# SPARSEHASH_INCLUDE_DIR - the SparseHash include directories

#if(SPARSEHASH_INCLUDE_DIR)
#set(SPARSEHASH_FIND_QUIETLY TRUE)
#endif(SPARSEHASH_INCLUDE_DIR)

find_path(SPARSEHASH_INCLUDE_DIR 
  sparsehash/internal/sparsehashtable.h
  PATHS ${SPARSEHASH_DIR}/include ${SPARSEHASH_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/external/
  NO_DEFAULT_PATH)

#find_path(SPARSEHASH_INCLUDE_DIR 
#  sparsehash/internal/sparsehashtable.h
#  PATHS ${CMAKE_CURRENT_SOURCE_DIR}/external/)

#if (NOT SPARSEHASH_INCLUDE_DIR)
#  find_path(SPARSEHASH_INCLUDE_DIR 
#  sparsehash/internal/sparsehashtable.h
#  PATHS ../external/)
#endif()

# handle the QUIETLY and REQUIRED arguments and set SPARSEHASH_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SparseHash DEFAULT_MSG SPARSEHASH_INCLUDE_DIR)

mark_as_advanced(SPARSEHASH_INCLUDE_DIR)

