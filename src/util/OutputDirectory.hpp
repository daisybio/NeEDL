//
// Created by juli on 22.06.22.
//

#ifndef GENEPISEEKER_OUTPUTDIRECTORY_HPP
#define GENEPISEEKER_OUTPUTDIRECTORY_HPP

#include <string>
#include <random>


namespace epi {

    class OutputDirectory {
    public:
        explicit OutputDirectory(std::string path, std::mt19937 random_device);

        std::string get_output_directory() const;

        std::ofstream get_ofstream(std::string name, std::string ending) const;
        std::string get_free_filepath(std::string name, std::string ending) const;
    private:
        std::string output_directory;
    };

} // epi


#ifdef HEADER_ONLY
#include "OutputDirectory.cpp"
#endif

#endif //GENEPISEEKER_OUTPUTDIRECTORY_HPP
