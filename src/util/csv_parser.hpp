/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of GenEpiSeeker.                                        *
 *                                                                             *
 *   GenEpiSeeker is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   GenEpiSeeker is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/



/*!
 * @file  csv_parser.hpp
 * @brief Declaration of epi::CSVParser.
 */

#ifndef SRC_UTIL_CSV_PARSER_HPP_
#define SRC_UTIL_CSV_PARSER_HPP_

#include "types.hpp"
#include "misc.hpp"

namespace epi {

/*!
 * @brief Parses CSV files.
 */
class CSVParser {

public:

	/*!
	 * @brief Constructor.
	 */
	CSVParser();

	/*!
	 * @brief Parse a CSV file.
	 * @param[in] filename Path to the CSV file.
	 * @param[in] sep Separator.
	 * @param[in] quote Quotation mark.
	 * @param[in] trim If @p true, leading and trailing white-spaces are trimmed.
	 */
	void parse(const std::string & filename, char sep = ',', char quote = '\"', bool trim = true, bool skip_empty_cells = true);

	/*!
	 * @brief Returns the number of rows in the CSV file.
	 * @return The number of rows.
	 */
	std::size_t num_rows() const;

	/*!
	 * @brief Returns the number of columns in the CSV file.
	 * @return The number of columns.
	 */
	std::size_t num_columns() const;

	/*!
	 * @brief Provides access to the parsed data.
	 * @param[in] row The ID of the row for which the data should be returned.
	 * @param[in] column The ID of the columns for which the data should be returned.
	 * @return The data contained in the cell (@p row, @p column).
	 */
	const std::string & cell(std::size_t row, std::size_t column) const;

private:

	std::size_t num_rows_;

	std::size_t num_columns_;

	std::vector<std::string> fields_;

	void parse_row_(const std::string & row, std::size_t row_id, char sep, char quote, bool trim, bool skip_empty_cells = true);

};

}

#ifdef HEADER_ONLY
#include "csv_parser.cpp"
#endif

#endif /* SRC_UTIL_CSV_PARSER_HPP_ */
