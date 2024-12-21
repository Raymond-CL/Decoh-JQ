cmake_minimum_required(VERSION 3.5)

project(
  decohJQ
  LANGUAGES Fortran
  DESCRIPTION "decoherent jet quenching"
  VERSION 0.5
)

enable_language(Fortran)

set(FVERSION "-std=legacy")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(fflag "-fno-strict-overflow")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} ${fflag}")

file(GLOB_RECURSE sources src/*.f90 src/*.f)

add_executable(eloss.exe ${sources} main.f90)