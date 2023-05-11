//
// Created by julian on 17.11.20.
//

#ifndef GENEPISEEKER_LOGGER_HPP
#define GENEPISEEKER_LOGGER_HPP

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include "types.hpp"
#include <boost/interprocess/sync/file_lock.hpp>

#define SECONDS_BETWEEN_PROGRESS_LOG 10

class Logger {
public:
    static void setConsoleLogging(bool consoleLogging);
    static void setFileLogging(bool fileLogging, std::string logDirectory);

    static void close();

    static bool consoleLoggingActive();

    static void logLine(std::string line);

    static void logCmdArgs(int argc, char** argv);

    static void logError(std::exception_ptr eptr);

    /*!
     * creates similar output as the linux /usr/bin/time -v command
     * only works under linux, will print empty message when used under other OS
     */
    static void logResourceUsage();

    /**
     * used to log a process (like local search) where updates should only be printed every 10 seconds
     * @param line
     */
    static void logProgress(std::string line);

    static void printHead();

    static std::ofstream* getOutputStream();


    template<class T>
    static std::string to_string(T val) {
        std::ostringstream out;
        out.precision(10);
        double val_d = double(val);
        if ((abs(val_d) < 0.0001 || abs(val_d) > 1000000.) && val_d != 0.) out << std::scientific;
        else out << std::fixed;
        out << val;
        return out.str();
    }

private:

    static void _printHead(std::ostream& out);
    static void _logLine(std::string line, std::ostream& out);
    static void _openFile(std::string logDirectory);
    static void _logCmdArgs();

    static bool _headPrinted;

    static bool _consoleLogging;
    static bool _fileLogging;
    static std::string _log_filepath;

    static std::ofstream *_log_filestream;

    static std::chrono::high_resolution_clock::time_point _lastProgressLog;

    static std::chrono::high_resolution_clock::time_point _resourceLoggingStartPoint;

    static std::vector<std::string> _cmdArgs;

};


#ifdef HEADER_ONLY
#include "Logger.cpp"
#endif

#endif //GENEPISEEKER_LOGGER_HPP
