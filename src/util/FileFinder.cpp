//
// Created by julian on 04.12.20.
//

#include "FileFinder.hpp"
#include <boost/filesystem.hpp>
#include "Logger.hpp"

using namespace boost::filesystem;

FileFinder::FileFinder(std::string error_info) {
    _error_info = error_info;
}

void FileFinder::add_starts_with(std::string prefix) {
    struct criterion crit {prefix, false,criterion_type::PREFIX};
    criterion_list.push_back(crit);
}

void FileFinder::add_ends_with(std::string suffix) {
    struct criterion crit {suffix, false,criterion_type::SUFFIX};
    criterion_list.push_back(crit);
}

void FileFinder::add_contains(std::string infix) {
    struct criterion crit {infix, false,criterion_type::INFIX};
    criterion_list.push_back(crit);
}

void FileFinder::add_not_starts_with(std::string prefix) {
    struct criterion crit {prefix, true,criterion_type::PREFIX};
    criterion_list.push_back(crit);
}

void FileFinder::add_not_ends_with(std::string suffix) {
    struct criterion crit {suffix, true,criterion_type::SUFFIX};
    criterion_list.push_back(crit);
}

void FileFinder::add_not_contains(std::string infix) {
    struct criterion crit {infix, true,criterion_type::INFIX};
    criterion_list.push_back(crit);
}

std::string FileFinder::get_file(std::string _directory) {
// iterate over all files in directory
    path dir(_directory);
    std::vector<std::string> selected_files; // To save the file names in a vector.

    if(is_directory(dir)) {
        for (auto it = directory_iterator(dir); it != directory_iterator(); ++it) {
            path path = it->path();
            std::string file = path.filename().string();
            if (is_regular_file(path)) {
                // test all criterions
                bool all_fulfilled = true;
                for (size_t i = 0; i < criterion_list.size(); i++) {
                    bool is_fulfilled = false;
                    switch(criterion_list[i].type) {
                        case PREFIX:
                            is_fulfilled = starts_with(file, criterion_list[i].content);
                            break;
                        case SUFFIX:
                            is_fulfilled = ends_with(file, criterion_list[i].content);
                            break;
                        case INFIX:
                            is_fulfilled = file.find(criterion_list[i].content) != std::string::npos;
                    }

                    if (is_fulfilled == criterion_list[i].isNot) {
                        // not fulfilled
                        all_fulfilled = false;
                        break;
                    }
                }

                if (all_fulfilled) {
                    selected_files.push_back(path.string());
                }
            }
        }
    }

    if (selected_files.size() == 0) {
        std::string curr_path = boost::filesystem::current_path().string();
        Logger::logLine("No match for file " + _error_info + " found");
        Logger::logLine("current working directory: " + curr_path);
        Logger::logLine("searched directory: " + _directory);
        throw epi::Error("File not found.");
    }

    if (selected_files.size() > 1) {
        std::string curr_path = boost::filesystem::current_path().string();
        Logger::logLine("Too many matches for file " + _error_info + " found");
        Logger::logLine("current working directory: " + curr_path);
        Logger::logLine("searched directory: " + _directory);
        throw epi::Error("File selection ambiguous.");
    }

    return selected_files[0];
}

bool FileFinder::starts_with(std::string str, std::string prefix) {
    if (str.length() < prefix.length()) return false;

    return str.substr(0, prefix.length()) == prefix;
}

bool FileFinder::ends_with(std::string str, std::string suffix) {
    if (str.length() < suffix.length()) return false;

    return str.substr(str.length() - suffix.length()) == suffix;
}
