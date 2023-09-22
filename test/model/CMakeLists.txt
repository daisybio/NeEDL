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


# Set executable directory.
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/test/model/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/test/model/bin)


# Add executables.
add_executable(compare_models src/compare_models.cpp)
target_link_libraries(compare_models boost_filesystem)
if(APPLE)
  target_link_libraries(compare_models omp)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_definitions(CMAKE_RELEASE)
  add_executable(NeEDL src/NeEDL.cpp src/util.cpp)
else()
  add_executable(NeEDL src/NeEDL.cpp src/util.cpp ${MAIN_SOURCES_NEW} ${QUANTUM_COMPUTING})
endif()

target_link_libraries(NeEDL sqlite3 boost_filesystem igraph ${Python_LIBRARIES})
if(APPLE)
  target_link_libraries(NeEDL omp)
endif()

add_executable(calculate_scores src/calculate_scores.cpp)
target_link_libraries(calculate_scores boost_filesystem igraph)
if(APPLE)
  target_link_libraries(calculate_scores omp)
endif()

add_executable(convert_to_binary src/convert_to_binary.cpp)
target_link_libraries(convert_to_binary boost_filesystem igraph)
if(APPLE)
  target_link_libraries(convert_to_binary omp)
endif()

add_executable(epiJSON src/epiJSON.cpp)
target_link_libraries(epiJSON boost_filesystem igraph)
if(APPLE)
  target_link_libraries(epiJSON omp)
endif()

# uWebSocekt stuff
file(GLOB uSocketObjects
        ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/*.o
)

add_executable(realtime_scores src/realtime_scores.cpp)
target_link_libraries(realtime_scores
        boost_filesystem
        igraph
        z
        m
        ${uSocketObjects}
        ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/lsquic/src/liblsquic/liblsquic.a
        ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/boringssl/build/ssl/libssl.a
        ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/boringssl/build/crypto/libcrypto.a
)
target_include_directories(realtime_scores
        PUBLIC ${PROJECT_SOURCE_DIR}/ext/uWebSockets/src
        PUBLIC ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/src
        PUBLIC ${PROJECT_SOURCE_DIR}/ext/uWebSockets/uSockets/boringssl/include
)

add_subdirectory(src)
