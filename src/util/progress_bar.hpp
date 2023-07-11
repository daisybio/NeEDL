/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of NeEDL.                                        *
 *                                                                             *
 *   NeEDL is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   NeEDL is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with NeEDL. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/

/*!
 * @file progress_bar.hpp
 * @brief ged::ProgressBar class declaration.
 */

#ifndef SRC_UTIL_PROGRESS_BAR_HPP_
#define SRC_UTIL_PROGRESS_BAR_HPP_

#include "types.hpp"

namespace epi {

/*!
 * @brief A progress bar class.
 */
class ProgressBar {

public:

	friend std::ostream & operator<<(std::ostream &, const ProgressBar &);

	/*!
	 * @brief Constructs a progress bar for given number of tasks.
	 * @param[in] num_tasks The number of tasks.
	 */
	ProgressBar(std::size_t num_tasks);

	/*!
	 * @brief Increments the number of solved tasks.
	 */
	void increment();

	/*!
	 * @brief Sets the number of solved tasks to 0.
	 */
	void reset();

	/*!
	 * @brief Sets the number of solved tasks to 0 and resets the number of tasks.
	 * @param[in] num_tasks The number of tasks.
	 */
	void reset(std::size_t num_tasks);

private:

	std::size_t num_solved_tasks_;

	std::size_t num_tasks_;

	std::chrono::high_resolution_clock::time_point start_time_;
};

/*!
 * @brief Streams the current progress in percent and the estimated remaining runtime.
 * @param[in,out] os Output stream.
 * @param[in] progress_bar The progress bar whose current progress should be streamed.
 * @return The output stream @p os.
 * @relates epi::ProgressBar
 */
std::ostream &
operator<<(std::ostream & os, const ProgressBar & progress_bar);

}

#ifdef HEADER_ONLY
#include "progress_bar.cpp"
#endif

#endif /* SRC_UTIL_PROGRESS_BAR_HPP_ */
