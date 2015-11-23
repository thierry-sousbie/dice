# Check availability of unordered maps
cmake_policy(PUSH)
cmake_minimum_required(VERSION 2.6.3)
cmake_policy(POP)

INCLUDE (CheckCXXSourceCompiles)
SET(CMAKE_CXX_FLAGS_BACKUP ${CMAKE_CXX_FLAGS})
UNSET(NATIVE_FLAGS_FOUND CACHE)
SET(NATIVE_CXX_FLAGS "-xHost")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BACKUP} ${NATIVE_CXX_FLAGS}")
check_cxx_source_compiles(
  "
#include <iostream>
  int main() 
   {
     double a=0;
     std::cout << a <<std::endl;
     return 0;
   }
"
  NATIVE_FLAGS_FOUND)


if (NOT NATIVE_FLAGS_FOUND) 
  
  UNSET(NATIVE_FLAGS_FOUND CACHE)  
  SET(NATIVE_CXX_FLAGS "-march=native")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BACKUP} ${NATIVE_CXX_FLAGS}")  
  check_cxx_source_compiles(
    "
#include <iostream>
  int main() 
   {
     double a=0;
     std::cout << a <<std::endl;
     return 0;
   }
  "
    NATIVE_FLAGS_FOUND)
endif()

if(NOT NATIVE_FLAGS_FOUND) 
  SET(NATIVE_CXX_FLAGS "")
endif()

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_BACKUP})

