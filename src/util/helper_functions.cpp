//
// Created by juli on 26.05.22.
//

#include "helper_functions.hpp"

#include <algorithm>
#include <boost/filesystem.hpp>
#include <regex>
#include "types.hpp"

namespace epi {
    std::string toUpperCase(std::string in) {
        std::transform(in.begin(), in.end(), in.begin(), ::toupper);
        return in;
    }

    std::string find_file_by_ending(std::string input_file, FileFinder constraints) {
        if (boost::filesystem::is_directory(input_file)) {
            input_file = constraints.get_file(input_file);
        } else if (!boost::filesystem::is_regular_file(input_file)) {
            throw epi::Error("File " + input_file + " not found.");
        }
        return input_file;
    }

    std::vector<std::string> string_split(const std::string &s, char delimiter, bool keep_empty) {
        std::vector<std::string> c;
        string_split(c, s, delimiter, keep_empty);
        return c;
    }

    std::string maskForJson(std::string input) {
        input = std::regex_replace(input, std::regex("\\\\"), "\\\\");
        return std::regex_replace(input, std::regex("\""), "\\\"");
    }

    std::string numberToBase32(unsigned long long input) {
        std::string result;

        int first_not_null = 0;

        const unsigned long long mask = 0x1f; // 5 x 1 in binary
        for (int i = 12; i >= 0; i--) {
            unsigned long long part = input >> std::min((i * 5), 59);
            if (part == 0) first_not_null++;
            part &= mask;

            if (part < 10) result += '0' + part;
            else result += 'A' + (part - 10);
        }

        result = result.substr(first_not_null);
        if (result.empty()) result = "0";
        return result;
    }

    std::string string_join(const std::vector<std::string> & vec, std::string delimiter) {
        std::stringstream ss;
        std::for_each(vec.begin(), vec.end(), [&ss, &delimiter] (const std::string &s) { ss << delimiter << s; });

        // ss.ignore(delimiter.size());
        auto str = ss.str();
        if (!str.empty()) str = str.substr(delimiter.size());

        return str;
    }


    unsigned long try_parse_number(const std::string & s, bool & successful) {
        char *pos;
        unsigned long v = strtoul(s.c_str(), &pos, 10);
        if (pos == 0) {
            successful = true;
            return v;
        }
        else {
            successful = false;
            return 0;
        }
    }


    std::unordered_map<std::string, char> getEscapedCharMap() {
        return {
            { "\\t", '\t' },
            { "\\n", '\n' },
            { "\\r", '\r' },
            { "-1", -1 }
            };
    }

    void run_plink_command(const std::string & ext_directory, std::vector<std::string> arguments) {
        arguments.insert(arguments.begin(), ext_directory + "plink/plink_linux_x86_64_20230116/plink");

        std::stringstream ss;
        for (auto arg : arguments) {
            boost::replace_all(arg, "\\", "\\\\");
            boost::replace_all(arg, "\"", "\\\"");
            ss << " \"" << arg << '"';
        }

        system(ss.str().c_str());
    }

    void remove_plink_temp_files(const std::string & path, std::vector<std::string> file_endings) {
        std::string log_msg = "Deleted temporary plink files with endings [";
        bool first_ending = true;
        for (auto &ending : file_endings) {
            if (boost::filesystem::exists(path + ending)) {
                boost::filesystem::remove(path + ending);
                if (first_ending) first_ending = false;
                else log_msg += ", ";
                log_msg += ending;
            }
        }
        log_msg += "] at '" + path + "'";

        Logger::logLine(log_msg);
    }



    double parseTimespanString(std::string s) {
        s.erase(std::remove(s.begin(), s.end(), ' '), s.end());

        double res_minutes = 0.;
        if (s.length() > 0) {
            char last = s.back();
            double factor = 1.;
            bool need_trunc = false;
            if (last == 's') {
                factor = 1.f/60.f;
                need_trunc = true;
            } else if (last == 'm') {
                factor = 1.;
                need_trunc = true;
            } else if (last == 'h') {
                factor = 60.;
                need_trunc = true;
            } else if (last == 'd') {
                factor = 1440.;
                need_trunc = true;
            }

            if (need_trunc) {
                s = s.substr(0, s.length() - 1);
            }

            res_minutes = std::stod(s);
            res_minutes *= factor;
        }
        return res_minutes;
    }

    void benjamini_hochberg_correction(std::vector<double> &pvalues) {
        if (pvalues.empty()) return;

        // convert to pvalue-index pairs
        std::vector<std::pair<double, size_t>> pval_index_pairs;
        for(size_t i = 0; i < pvalues.size(); i++) {
            pval_index_pairs.push_back({ pvalues[i], i });
        }

        // sort pairs in ascending order by p value
        std::sort(pval_index_pairs.begin(), pval_index_pairs.end(), [] (std::pair<double, size_t> a, std::pair<double, size_t> b) { return a.first < b.first; });

        // correct p-value
        for (size_t i = 0; i < pval_index_pairs.size(); i++) {
            pval_index_pairs[i].first *= double(pval_index_pairs.size()) / double(i + 1);
        }

        // check that lower p-value cannot lead to higher corrected one comparing the others
        double min = pval_index_pairs.back().first;
        for (long i = pval_index_pairs.size() - 1; i >= 0; i--) {
            if (pval_index_pairs[i].first > min) {
                pval_index_pairs[i].first = min;
            } else {
                min = pval_index_pairs[i].first;
            }
        }

        // sort back adjusted p-values
        for (auto x : pval_index_pairs) {
            pvalues[x.second] = x.first;
        }
    }

} // epi