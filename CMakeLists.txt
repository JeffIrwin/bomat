
cmake_minimum_required(VERSION 3.7)

#set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_COMPILER "ifort")

set(CMAKE_C_COMPILER "icc")
set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_EXE_LINKER "ifort")
set(CMAKE_C_LINK_EXECUTABLE "ifort")
set(CMAKE_LINKER "ifort")

set(PROJECT bomat)

project(${PROJECT} Fortran CXX)
#project(${PROJECT} Fortran)

# TODO:  build LAPACK if necessary via make

# TODO: remove uninit flag.  It's for a gfortran bug
set(CMAKE_Fortran_FLAGS "-W:all -W:extra -qopenmp -qmkl -fPIC")
#set(CMAKE_Fortran_FLAGS "-Wall -Wextra -fopenmp -Wno-tabs -Wno-maybe-uninitialized")

set(CMAKE_Fortran_FLAGS_DEBUG "-pedantic -fbounds-check -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

if(APPLE)
	set(CMAKE_EXE_LINKER_FLAGS "-openmp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -qmkl")
endif()

message("")
message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
message("")

set(SUBMODULE_DIR submodules)

#set(PNMIO pnmio)
#set(JSONFORTRAN json-fortran)
set(COLORMAPPER colormapper)

set(LAPACK_DIR ${SUBMODULE_DIR}/lapack)

#add_subdirectory(${SUBMODULE_DIR}/${PNMIO})
#add_subdirectory(${SUBMODULE_DIR}/${JSONFORTRAN})
add_subdirectory(${SUBMODULE_DIR}/${COLORMAPPER})

add_executable(${PROJECT} main.f90 colormapper_wrapper.cpp)

set_property(TARGET ${PROJECT} PROPERTY LINKER_LANGUAGE Fortran)

#include_directories(
##	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${PNMIO}
##	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${JSONFORTRAN}/include
#	)

target_link_libraries(${PROJECT}
	#${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/liblapack.a
	#${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/librefblas.a
	${COLORMAPPER}
#	#${PNMIO}
#	#jsonfortran-static
	)

