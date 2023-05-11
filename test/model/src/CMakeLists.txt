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
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/test/model/bin)

add_executable(create_PoCos_input create_PoCos_input.cpp)
if(APPLE)
    target_link_libraries(create_PoCos_input omp)
endif()
