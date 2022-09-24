cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_COMPILER "icc")
set(CMAKE_CXX_COMPILER "icpc")
project(FE_1D)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fma -AVX2 -std=c++14 -lmkl_rt -fopenmp")

include_directories(${CMAKE_SOURCE_DIR}/include/)
include_directories(${CMAKE_SOURCE_DIR}/include/UMFPACK/)

add_subdirectory(src)

# Set 1 for DEBUG mode and 0 for normal Execution
set(DEBUG_FLAG "0" CACHE STRING "Enter the Debug Flag" )

#set 1 for problem type 1 ( FEM )  2 for others
set(PROB_TYPE "1" CACHE STRING "Enter the Debug Flag" )




if(${PROB_TYPE} EQUAL "1")
    add_executable(prob1 MainPrograms/Prob_1.cpp)
    target_link_libraries(prob1 FE_1D ${CMAKE_SOURCE_DIR}/lib/umfpack.a)
endif()

if(${PROB_TYPE} EQUAL "2")
    add_executable(prob2 MainPrograms/Prob_2.cpp)
    target_link_libraries(prob2 FE_1D ${CMAKE_SOURCE_DIR}/lib/umfpack.a)
endif()

set(CMAKE_EXE_LINKER_FLAGS "-O2 -fma -AVX2 -ldl -lpthread")

if(${DEBUG_FLAG} EQUAL "1")
  add_definitions( -D_DEBUG  )
  message("-------------- Debug Flag Set  ----------" )
endif() 





