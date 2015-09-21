# Check availability of unordered maps
cmake_policy(PUSH)
cmake_minimum_required(VERSION 2.6.3)
cmake_policy(POP)

INCLUDE (CheckCXXSourceCompiles)
UNSET(CPP11_FOUND)
SET(CMAKE_CXX_FLAGS_BACKUP ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #-Wcast-qual -Weffc++
check_cxx_source_compiles(
  "
  #include <unordered_map>
  int main() 
   {
     std::unordered_map<int, int> m;
     return 0;
   }
  "
CPP11_FOUND)
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_BACKUP})

