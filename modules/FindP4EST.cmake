find_path(P4EST_INCLUDE_DIR 
  p4est.h
  PATHS ${P4EST_DIR}/include ${P4EST_DIR} 
  NO_DEFAULT_PATH)

find_path(P4EST_INCLUDE_DIR p4est.h)

find_library(P4EST_LIB_DIR p4est
  PATHS ${P4EST_DIR}/lib ${P4EST_DIR} 
  NO_DEFAULT_PATH)

find_library(P4EST_LIB_DIR p4est)

if (P4EST_LIB_DIR AND P4EST_INCLUDE_DIR)
  SET(P4EST_FOUND true)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(p4est DEFAULT_MSG P4EST_INCLUDE_DIR P4EST_LIB_DIR)

mark_as_advanced(P4EST_INCLUDE_DIR P4EST_LIB_DIR)
