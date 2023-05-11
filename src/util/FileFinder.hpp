//
// Created by julian on 04.12.20.
//

#ifndef GENEPISEEKER_FILEFINDER_H
#define GENEPISEEKER_FILEFINDER_H

#include <string>
#include <fstream>
#include <vector>

/**
 * Helper class to select a distinct file in a directory based on provided criterions.
 * If more than one file exists in directory with the given criterions an error is thrown.
 */
class FileFinder {

public:
    FileFinder() = default;
    FileFinder(std::string error_info);

    void add_starts_with(std::string prefix);
    void add_ends_with(std::string suffix);
    void add_contains(std::string infix);

    void add_not_starts_with(std::string prefix);
    void add_not_ends_with(std::string suffix);
    void add_not_contains(std::string infix);

    std::string get_file(std::string directory);

private:
    std::string _error_info;

    enum criterion_type {
        PREFIX, SUFFIX, INFIX
    };

    struct criterion{
        std::string content;
        bool isNot = false;
        criterion_type type;
    };

    std::vector<struct criterion> criterion_list;

    bool starts_with(std::string str, std::string prefix);
    bool ends_with(std::string str, std::string suffix);
};

#ifdef HEADER_ONLY
#include "FileFinder.cpp"
#endif

#endif //GENEPISEEKER_FILEFINDER_H