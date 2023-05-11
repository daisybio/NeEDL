//
// Created by juli on 22.06.22.
//

#include "OutputDirectory.hpp"
#include <chrono>
#include <iomanip>
#include <utility>
#include "types.hpp"

namespace epi {
    OutputDirectory::OutputDirectory(std::string path, std::mt19937 random_device) {
        if (path.back() != '/' && path.back() != '\\') path += '/';

        // get time string as prefix for folder
        auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::tm *ltime = std::localtime(&now);
        std::stringstream ss;
        ss << std::put_time(ltime, "%y_%m_%d-%H_%M_%S");
        std::string time_str = ss.str();

        std::string dir_path;
        std::uniform_int_distribution<unsigned long long> dist(0, std::numeric_limits<unsigned long long>::max());
        do {
            unsigned long long random_val = dist(random_device);
            dir_path = path + time_str + "-" + std::to_string(random_val);
        } while (boost::filesystem::exists(dir_path));

        boost::filesystem::create_directory(dir_path);

        output_directory = dir_path + '/';
    }

    std::string OutputDirectory::get_output_directory() const {
        return output_directory;
    }

    std::ofstream OutputDirectory::get_ofstream(std::string name, std::string ending) const {
        std::string path = get_free_filepath(std::move(name), std::move(ending));
        return std::ofstream(path);
    }

    std::string OutputDirectory::get_free_filepath(std::string name, std::string ending) const {
        std::string path = output_directory + name + ending;
        unsigned long index = 2;
        while (boost::filesystem::exists(path)) {
            path = output_directory + name;
            path += '_';
            path += std::to_string(index);
            path += ending;
            index ++;
        }

        return path;
    }
} // epi