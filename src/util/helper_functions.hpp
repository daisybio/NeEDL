//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_HELPER_FUNCTIONS_HPP
#define GENEPISEEKER_HELPER_FUNCTIONS_HPP

#include <string>
#include <vector>
#include <algorithm>
#include "FileFinder.hpp"
#include "types.hpp"

namespace epi {

    std::string toUpperCase(std::string in);

    template<class T>
    bool includes(const std::vector<T>& haystack, T needle) {
        for(auto &item : haystack) {
            if (item == needle) return true;
        }
        return false;
    }

    /*
    template<class T>
    void sort_omp(std::vector<T>& data) {
        size_t num_blocks = omp_get_max_threads();
        if (data.size() < 10000 || num_blocks == 1) {
            std::sort(data.begin(), data.end());
        } else {
            // sort blocks individually
#pragma omp parallel for default(none) shared(data, num_blocks) schedule(static)
            for (size_t i = 0; i < num_blocks; i++) {
                size_t start = data.size() * i / num_blocks;
                size_t end = data.size() * (i + 1) / num_blocks;
                std::sort(data.begin() + start, data.begin() + end);
            }

            // merge blocks
            for (size_t i = 0; i < num_blocks; ++i) {
                size_t middle = data.size() * i / num_blocks;
                size_t end = data.size() * (i + 1) / num_blocks;
                std::inplace_merge(data.begin(), data.begin() + middle, data.begin() + end);
            }
        }
    }
     */


    std::string find_file_by_ending(std::string path, FileFinder constraints);

    // performance optimization (faster than boost method)
    template<typename Container>
    void string_split(Container &results, const typename Container::value_type & s, typename Container::value_type::value_type delimiter, bool keep_empty = true) {
        results.clear();
        std::istringstream ss(s);
        typename Container::value_type field;
        while (!ss.eof()) {
            getline(ss, field, delimiter);
            if (keep_empty || !field.empty())
                results.push_back(field);
        }
    }

    // wrapper of the above function which creates a new container
    std::vector<std::string> string_split(const std::string & s, char delimiter, bool keep_empty = true);


    std::string maskForJson(std::string input);

    std::string numberToBase32(unsigned long long input);

    std::string string_join(const std::vector<std::string> & vec, std::string delimiter);

    template<class T>
    std::string string_join(const std::vector<T> & vec, std::string delimiter, std::string (*converter) (T)) {
        std::stringstream ss;
        std::for_each(vec.begin(), vec.end(), [&ss, &delimiter, &converter] (const T &s) { ss << delimiter << converter(s); });

        // ss.ignore(delimiter.size());
        auto str = ss.str();
        if (!str.empty()) str = str.substr(delimiter.size());

        return str;
    }

    void benjamini_hochberg_correction(std::vector<double> &pvalues);


    std::unordered_map<std::string, char> getEscapedCharMap ();

    unsigned long try_parse_number(const std::string & s, bool & successful);


    void run_plink_command(const std::string & ext_directory, std::vector<std::string> arguments);
    void remove_plink_temp_files(const std::string & path, std::vector<std::string> file_endings);


    double parseTimespanString(std::string s);

    template <class T>
    void expand_multi_param(std::string name, std::vector<T> & params, T default_value, size_t size, bool is_mandatory) {
        if (size == 0) {
            params.clear();
            return;
        }

        if (params.empty()) {
            if (is_mandatory) throw epi::Error("Parameter " + name + " is mandatory but no values were provided.");
            else params.push_back(default_value);
        }

        if (params.size() == 1) {
            // expand to size
            for (size_t i = 1; i < size; i++) {
                params.push_back(params[0]);
            }
        } else if (params.size() != size) {
            // not correct number of params
            throw epi::Error("Parameter " + name + " was given " + std::to_string(params.size()) + " times, but either 1 or " + std::to_string(size) + " were expected.");
        }
    }
} // epi

#ifdef HEADER_ONLY
#include "helper_functions.cpp"
#endif

#endif //GENEPISEEKER_HELPER_FUNCTIONS_HPP
