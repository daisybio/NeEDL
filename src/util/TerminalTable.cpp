//
// Created by julian on 18.12.20.
//

#include "TerminalTable.hpp"

#include <iostream>


TerminalTable::TerminalTable(size_t trunc_max_width, size_t trunc_max_rows, size_t spacing, size_t min_width) {
    this->trunc_max_width = trunc_max_width;
    this->trunc_max_rows = trunc_max_rows;
    this->spacing = spacing;
    this->min_width = min_width;
}

void TerminalTable::add_row(std::vector<std::string> cells) {
    table.push_back(cells);
    if (col_widths.size() < cells.size()) {
        size_t missing = cells.size() - col_widths.size();
        for (size_t i = 0; i < missing; i++) {
            col_widths.push_back(min_width + spacing);
        }
    }

    for(size_t i = 0; i < cells.size(); i++) {
        col_widths[i] = std::max(cells[i].size() + spacing, col_widths[i]);
    }
}

void TerminalTable::print() {
    size_t included_columns = col_widths.size();
    if (trunc_max_width != 0) {
        // width limit is set --> check how many columns can be displayed before the limit is reached
        size_t total_width = 0;
        for (included_columns = 0; included_columns < col_widths.size(); ++included_columns) {
            total_width += col_widths[included_columns];
            if (total_width > trunc_max_width) break;
        }
    }

    size_t included_rows = table.size();
    if (trunc_max_rows != 0) {
        included_rows = std::min(trunc_max_rows + 1, table.size());
    }

    for (size_t i = 0; i < included_rows; i++) {
        for (size_t j = 0; j < std::min(table[i].size(), included_columns); j++) {

            std::cout << table[i][j];
            int spaces = col_widths[j] - table[i][j].size();
            if (spaces > 0) {
                std::cout << std::string(spaces, ' ');
            }
        }

        if (included_columns < col_widths.size()) {
            if (i == 0) {
                std::cout << " ... (+ " << (col_widths.size() - included_columns) << " cols)";
            } else {
                std::cout << " ...";
            }
        }

        std::cout << std::endl;
        if (i == 0) {
            size_t total_width = 0;
            for (size_t j = 0; j < included_columns; j++) {
                total_width += col_widths[j];
            }
            if (included_columns < col_widths.size()) {
                // continue line under truncate msg
                total_width += 14 + std::to_string(col_widths.size() - included_columns).size();
            }
            std::cout << std::string(total_width, '=') << std::endl;
        }
    }

    if (included_rows < table.size()) {
        std::cout << "... (+ " << (table.size() - included_rows) << " rows)" << std::endl;
    }
}
