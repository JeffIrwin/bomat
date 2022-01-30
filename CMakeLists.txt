
cmake_minimum_required(VERSION 3.7)

set(PROJECT bomat)

project(${PROJECT} Fortran CXX)

# TODO:  build LAPACK if necessary via make

# TODO: remove uninit flag.  It's for a gfortran bug
set(CMAKE_Fortran_FLAGS "-Wall -Wextra -fopenmp -Wno-tabs -Wno-maybe-uninitialized")

set(CMAKE_Fortran_FLAGS_DEBUG "-pedantic -fbounds-check -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

if(APPLE)
	set(CMAKE_EXE_LINKER_FLAGS "-openmp -Wl,-no_compact_unwind")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")
endif()

set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(SUBMODULE_DIR submodules)

#set(PNMIO pnmio)
#set(JSONFORTRAN json-fortran)
set(COLORMAPPER colormapper)

set(LAPACK_DIR ${SUBMODULE_DIR}/lapack)

#add_subdirectory(${SUBMODULE_DIR}/${PNMIO})
#add_subdirectory(${SUBMODULE_DIR}/${JSONFORTRAN})
add_subdirectory(${SUBMODULE_DIR}/${COLORMAPPER})

add_executable(${PROJECT}
	bomat.f90
	colormapper_wrapper.cpp
	)

#include_directories(
##	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${PNMIO}
##	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${JSONFORTRAN}/include
#	)

target_link_libraries(${PROJECT}
	${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/liblapack.a
	${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/librefblas.a
	${COLORMAPPER}
#	#${PNMIO}
#	#jsonfortran-static
	)

