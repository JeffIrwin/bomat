
cmake_minimum_required(VERSION 3.7)

if(DEFINED ENV{BOMAT_INTEL})

	set(CMAKE_Fortran_COMPILER "ifort")
	set(CMAKE_C_COMPILER "icc")
	set(CMAKE_CXX_COMPILER "icpc")

else()

	set(CMAKE_Fortran_COMPILER "gfortran")

endif()

set(PROJECT bomat)

project(${PROJECT} Fortran CXX)

# TODO:  build LAPACK if necessary via make for gfortran (no need for ifort)

if(DEFINED ENV{BOMAT_INTEL})
	set(CMAKE_Fortran_FLAGS "-W1 -qopenmp -qmkl -fPIC")
else()
	set(CMAKE_Fortran_FLAGS "-Wall -Wextra -fopenmp -Wno-tabs -Wno-maybe-uninitialized")
endif()

message("")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message("")

set(CMAKE_Fortran_FLAGS_DEBUG "-pedantic -fbounds-check -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

if(APPLE)

	set(CMAKE_EXE_LINKER_FLAGS "-openmp")
	#set(CMAKE_EXE_LINKER_FLAGS "-openmp -Wl,-no_compact_unwind")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	#set(CMAKE_C_COMPILER "gcc")
	#set(CMAKE_CXX_COMPILER "gcc")
	#set(CMAKE_CXX_COMPILER "g++")

else()

	if(DEFINED ENV{BOMAT_INTEL})
		set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -qmkl")
	else()
		set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")
	endif()

endif()

message("")
message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
message("")

set(SUBMODULE_DIR submodules)

#set(PNMIO pnmio)
set(JSONFORTRAN json-fortran)
set(COLORMAPPER colormapper)

set(LAPACK_DIR ${SUBMODULE_DIR}/lapack)

#add_subdirectory(${SUBMODULE_DIR}/${PNMIO})
add_subdirectory(${SUBMODULE_DIR}/${JSONFORTRAN})
add_subdirectory(${SUBMODULE_DIR}/${COLORMAPPER})

add_executable(${PROJECT}
	bomat.f90
	colormapper_wrapper.cpp
	)

# This tells cmake to use the fortran compiler as the linker.  See:
#
#     https://stackoverflow.com/questions/40767212/use-fortran-linker-in-cmake-mixed-language-project
#
set_property(TARGET ${PROJECT} PROPERTY LINKER_LANGUAGE Fortran)

include_directories(
##	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${PNMIO}
	${CMAKE_CURRENT_BINARY_DIR}/${SUBMODULE_DIR}/${JSONFORTRAN}/include
	)

set(LINK_LIBS
	${COLORMAPPER}
#	#${PNMIO}
	jsonfortran-static
	)

if(NOT DEFINED ENV{BOMAT_INTEL})

	list(APPEND LINK_LIBS
		${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/liblapack.a
		${CMAKE_SOURCE_DIR}/${LAPACK_DIR}/librefblas.a
		)

endif()

target_link_libraries(${PROJECT} ${LINK_LIBS})

message("")
message("LINK_LIBS = ${LINK_LIBS}")
message("")

