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

#include "csv_parser.hpp"

/*!
 * @file  csv_parser.ipp
 * @brief Definition of epi::CSVParser.
 */

#ifndef SRC_UTIL_CSV_PARSER_IPP_
#define SRC_UTIL_CSV_PARSER_IPP_

namespace epi {

CSVParser::
CSVParser():
num_rows_{undefined_uint()},
num_columns_{undefined_uint()},
fields_() {}

void
CSVParser::
parse(const std::string & filename, char sep, char quote, bool trim, bool skip_empty_cells) {
	num_rows_ = 0;
	num_columns_ = undefined_uint();
	fields_.clear();
	std::ifstream csv_file(filename);
	if (not csv_file.is_open()) {
		throw Error("The file " + filename + " cannot be opened.");
	}
	std::string row;
	std::size_t row_id{0};
	while(std::getline(csv_file, row)) {
		parse_row_(row, row_id++, sep, quote, trim, skip_empty_cells);
	}
	csv_file.close();
}

std::size_t
CSVParser::
num_rows() const {
	return num_rows_;
}

std::size_t
CSVParser::
num_columns() const {
	return num_columns_;
}

const std::string &
CSVParser::
cell(std::size_t row, std::size_t column) const {
	if (row >= num_rows_) {
		throw Error(std::string("Non-existing row ") + std::to_string(row) + ". Must be between 0 and " + std::to_string(num_rows_ - 1) + ".");
	}
	if (column >= num_columns_) {
		throw Error(std::string("Non-existing column ") + std::to_string(column) + ". Must be between 0 and " + std::to_string(num_columns_ - 1) + ".");
	}
	return fields_.at(row * num_columns_ + column);
}

void
CSVParser::
parse_row_(const std::string & row, std::size_t row_id, char sep, char quote, bool trim, bool skip_empty_cells) {
	std::size_t num_fields_in_row{misc::tokenize(row, sep, quote, trim, fields_, false, skip_empty_cells)};
	if (num_columns_ == undefined_uint()) {
		num_columns_ = num_fields_in_row;
	}
	if (num_columns_ != num_fields_in_row) {
		throw Error("The rows " + std::to_string(row_id - 1) + " and " + std::to_string(row_id) + " have different numbers of columns.");
	}
	num_rows_++;
}


}


#endif /* SRC_UTIL_CSV_PARSER_IPP_ */
