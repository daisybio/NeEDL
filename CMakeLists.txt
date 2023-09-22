###############################################################################
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of NeEDL.                                        #
#                                                                             #
#   NeEDL is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   NeEDL is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with NeEDL. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
###############################################################################

cmake_minimum_required(VERSION 3.8)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(NeEDL)

include(CMakePrintHelpers)

# Determine build type.
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Building NeEDL with build type 'Release', as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release or Debug build type")
else()
  message(STATUS "Building NeEDL with build type '${CMAKE_BUILD_TYPE}'.")
endif()

set(USE_INCLUDED_BOOST "1" CACHE BOOL "Use the boost version shipped with this repository")
cmake_print_variables(USE_INCLUDED_BOOST)
if (USE_INCLUDED_BOOST)
  message("-- boost version: vendored")
  set(BOOST_INCLUDE ${CMAKE_SOURCE_DIR}/ext/boost_1_71_0)
  set(BOOST_LINK_DIR ${CMAKE_SOURCE_DIR}/ext/boost_1_71_0/stage/lib)
else()
  message("-- boost version: system")
  set(Boost_USE_STATIC_LIBS OFF) 
  set(Boost_USE_MULTITHREADED ON)  
  set(Boost_USE_STATIC_RUNTIME OFF) 
  find_package(Boost 1.71.0 REQUIRED COMPONENTS filesystem graph)

  set(BOOST_INCLUDE ${Boost_INCLUDE_DIRS})
endif()

# Find doxygen.
#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
#  add_custom_target(docs ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Generating documentation with Doxygen." VERBATIM)
#endif()

# needed for the quantum-computing part
find_package(Python COMPONENTS Interpreter Development)
find_package(SQLite3 REQUIRED)

# debug stuff
cmake_print_variables(CMAKE_C_COMPILER)
cmake_print_variables(CMAKE_C_COMPILER_VERSION)
cmake_print_variables(CMAKE_CXX_COMPILER)
cmake_print_variables(CMAKE_CXX_COMPILER_VERSION)
cmake_print_variables(SQLite3_INCLUDE_DIRS)
cmake_print_variables(SQLite3_LIBRARIES)

# Set up the compiler.
set(CMAKE_CXX_FLAGS "-std=c++20 -Wall -Wextra -pedantic -Wredundant-decls -Wdisabled-optimization -Wctor-dtor-privacy -Wnon-virtual-dtor -Werror=old-style-cast -Wsign-promo -Werror=return-type -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(APPLE)
  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-instr-generate -fcoverage-mapping")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif()
set(CMAKE_MACOSX_RPATH ON)

# Add include and link directories.
set(EIGEN_ROOT ${CMAKE_SOURCE_DIR}/ext/eigen_3_3_7)
set(CLI11_ROOT ${CMAKE_SOURCE_DIR}/ext/cli11_1_9_0)
set(CATCH_ROOT ${CMAKE_SOURCE_DIR}/ext/catch_2_13_9)
set(IGRAPH_ROOT ${CMAKE_SOURCE_DIR}/ext/igraph_0.9.8)
set(PYBIND11_HOME ${CMAKE_SOURCE_DIR}/ext/pybind11)

include_directories(SYSTEM ${SQLite3_INCLUDE_DIRS} ${BOOST_INCLUDE} ${EIGEN_ROOT} ${CLI11_ROOT} ${CATCH_ROOT} ${IGRAPH_ROOT}/include ${IGRAPH_ROOT}/build/include ${PYBIND11_HOME}/include ${SQLite3_INCLUDE_DIRS} ${Python_INCLUDE_DIRS} /usr/lib/llvm-10/include)
link_directories(${BOOST_LINK_DIR} /usr/lib/llvm-10/lib ${IGRAPH_ROOT}/build/src ${Python_LIBRARY_DIRS} ${SQLite3_LIBRARIES})
if(APPLE)
  include_directories(SYSTEM ${OMP_ROOT}/include)
  link_directories(${OMP_ROOT}/lib)
endif()

# add main sources
file(GLOB_RECURSE MAIN_SOURCES_NEW src/*.h src/*.cpp src/*.hpp)
file(GLOB_RECURSE QUANTUM_COMPUTING quepistasis/header/*.h quepistasis/src/*.cpp quepistasis/header/*.hpp)

# Add subdirectories.
add_subdirectory(ext)
add_subdirectory(test)
add_subdirectory(src)
