# Check advanced features
cmake_policy(PUSH)
cmake_minimum_required(VERSION 2.6.3)
cmake_policy(POP)

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(cmessage)
UNSET(VALID_COMPILE_FLAGS_FOUND CACHE)
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
VALID_COMPILE_FLAGS_FOUND)
if(NOT VALID_COMPILE_FLAGS_FOUND)
cmessage(STATUS_RED " ERROR: invalid compilation flags: ${CMAKE_CXX_FLAGS}.")
message(FATAL_ERROR "Cannot continue ...")
endif()

SET(CMAKE_CXX_FLAGS_BACKUP ${CMAKE_CXX_FLAGS})
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #-Wcast-qual -Weffc++
UNSET(NULLPTR_FOUND CACHE)
check_cxx_source_compiles(
  "
  int main() 
   {
     double *a=nullptr;
     return 0;
   }
  "
NULLPTR_FOUND)

UNSET(OMP_ATOMIC_CAPTURE_FOUND CACHE)
check_cxx_source_compiles(
  "
  int main() 
   {
    int i=0;
    int j=0;
    #pragma omp parallel
    for (int i=0;i<10;)
    {
    #pragma omp atomic capture
     { 
      j=i;
      i+=1;
     }
    }
    return 0;
   }
  "
OMP_ATOMIC_CAPTURE_FOUND)

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_BACKUP})

