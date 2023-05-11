//
// Created by julian on 18.12.20.
//

#ifndef GENEPISEEKER_TERMINALTABLE_H
#define GENEPISEEKER_TERMINALTABLE_H

#include <string>
#include <vector>

class TerminalTable {
public:
    TerminalTable(size_t trunc_max_width = 0, size_t trunc_max_rows = 0, size_t spacing = 3, size_t min_width = 5);

    void add_row(std::vector<std::string> cells);
    void print();
private:
    std::vector<std::vector<std::string>> table;
    std::vector<size_t> col_widths;
    size_t trunc_max_width;
    size_t trunc_max_rows;
    size_t min_width;
    size_t spacing;
};


#ifdef HEADER_ONLY
#include "TerminalTable.cpp"
#endif

#endif //GENEPISEEKER_TERMINALTABLE_H
