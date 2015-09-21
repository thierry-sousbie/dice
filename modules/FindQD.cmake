

# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# FIXME: How do I find the version of QD that I want to use?
# What versions are available?

# NOTE: QD prefix is understood to be the path to the root of the QD
# installation library.
set(QD_PREFIX "" CACHE PATH "The path to the prefix of a QD installation")
set(QD_DIR "" CACHE PATH "The path to the root of a QD installation")

find_path(QD_INCLUDE_DIR qd/qd_real.h
PATHS ${QD_DIR} ${QD_DIR}/include ${QD_PREFIX} ${QD_PREFIX}/include /usr/include /usr/local/include)

find_library(QD_LIBRARY NAMES qd
PATHS ${QD_DIR} ${QD_DIR}/lib ${QD_PREFIX} ${QD_PREFIX}/lib /usr/lib /usr/local/lib)

if(QD_INCLUDE_DIR AND QD_LIBRARY)
get_filename_component(QD_LIBRARY_DIR ${QD_LIBRARY} PATH)
set(QD_FOUND TRUE)
endif()

if(QD_FOUND)
if(NOT QD_FIND_QUIETLY)
MESSAGE(STATUS "Found QD: ${QD_LIBRARY}")
endif()
elseif(QD_FOUND)
if(QD_FIND_REQUIRED)
message(FATAL_ERROR "Could not find QD")
endif()
endif()
