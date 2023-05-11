###############################################################################
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of GenEpiSeeker.                                        #
#                                                                             #
#   GenEpiSeeker is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GenEpiSeeker is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
###############################################################################

# Set executable directory.
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/test/unit/bin)

# Custom target for all unit tests.
add_custom_target(unit DEPENDS instance variance_model penetrance_model regression_model bayesian_model)

# Add executables.
add_executable(instance src/instance.cpp src/catch_main.cpp)
add_executable(variance_model src/variance_model.cpp src/catch_main.cpp)
add_executable(penetrance_model src/penetrance_model.cpp src/catch_main.cpp)
add_executable(regression_model src/regression_model.cpp src/catch_main.cpp)
if(APPLE)
target_link_libraries(regression_model omp)
endif()
add_executable(bayesian_model src/bayesian_model.cpp src/catch_main.cpp)