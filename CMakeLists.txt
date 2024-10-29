cmake_minimum_required(VERSION 3.5)

project(
  eloss
  LANGUAGES Fortran
  VERSION 0.1
)

enable_language(Fortran)

set(FVERSION "-std=f95")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(fflag "-fno-strict-overflow")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${fflag}")

file(GLOB_RECURSE sources src/*.f90 src/*.f)

add_executable(eloss ${sources} main.f90)