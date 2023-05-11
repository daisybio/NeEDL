//
// Created by julian on 17.11.20.
//

#include "Logger.hpp"

#ifdef __linux__
#include <sys/time.h>
#include <sys/resource.h>
#endif

bool Logger::_headPrinted = false;
bool Logger::_consoleLogging = true;
bool Logger::_fileLogging = false;
std::ofstream *Logger::_log_filestream = nullptr;
std::string Logger::_log_filepath { "" };
std::chrono::high_resolution_clock::time_point Logger::_resourceLoggingStartPoint = std::chrono::high_resolution_clock::now();
std::chrono::high_resolution_clock::time_point Logger::_lastProgressLog = std::chrono::high_resolution_clock::now();
std::vector<std::string> Logger::_cmdArgs;

void Logger::setConsoleLogging(bool consoleLogging) {
    _consoleLogging = consoleLogging;
}

void Logger::_openFile(std::string logDirectory) {
    if (_log_filestream != nullptr) delete _log_filestream;
    _log_filepath = logDirectory + "run.log";

    _log_filestream = new std::ofstream(_log_filepath);
}

void Logger::setFileLogging(bool fileLogging, std::string logDirectory) {
    _fileLogging = fileLogging;
    if (fileLogging) _openFile(logDirectory);
}

void Logger::close() {
    if (_log_filestream != nullptr) {
        _log_filestream->close();
    }
}

void Logger::logLine(std::string line) {
    if (!_headPrinted) {
        printHead();
    }

    if (_consoleLogging) _logLine(line, std::cout);
    if (_fileLogging) {
        _logLine(line, *_log_filestream);
    }
}

void Logger::_logLine(std::string line, std::ostream& out) {
    auto now_exact = std::chrono::system_clock::now();
    auto since_epoch = now_exact.time_since_epoch();
    auto micros = std::chrono::duration_cast<std::chrono::microseconds>(since_epoch).count();
    auto now = std::chrono::system_clock::to_time_t(now_exact);
    std::tm *ltime = std::localtime(&now);
    long micros_diff = micros - (now * 1000000L);
    std::string micros_str = std::to_string(micros_diff);
    out << "[" << std::put_time(ltime, "%y/%m/%d-%H:%M:%S") << '.';
    for (size_t i = micros_str.size(); i < 6; i++) out << '0';
    out << micros_diff << "]  " << line << std::endl;
}

void Logger::printHead() {
    if (!_headPrinted) {
        if (_consoleLogging) {
            _printHead(std::cout);
        }
        if (_fileLogging) {
            _printHead(*_log_filestream);
        }
        _headPrinted = true;

        if(!_cmdArgs.empty()) _logCmdArgs();
        logLine("Log file: " + (_fileLogging ? _log_filepath : "---"));
    }
}

void Logger::_printHead(std::ostream& out) {
    out << "***************************************************************\n";
    out << "                            NeEDL                              \n";
    out << "  Scalable network-based epistasis detection via local search  \n";
    out << "                          version 1.0                          \n";
    out << "***************************************************************\n";
}

void Logger::logProgress(std::string line) {
    std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
    if (std::chrono::duration_cast<std::chrono::seconds>(now - _lastProgressLog).count() >= SECONDS_BETWEEN_PROGRESS_LOG) {
        _lastProgressLog = now;
        logLine(line);
    }
}

bool Logger::consoleLoggingActive() {
    return _consoleLogging;
}



void Logger::logResourceUsage() {
#ifdef __linux__
    struct rusage usage_info;
    int state = getrusage(RUSAGE_SELF, &usage_info);

    if (state != 0) {
        logLine("An error occurred when retrieving resource usage information.");
        return;
    }

    // print usage information to logfile
    logLine("RUN STATS: ");

    long long user_time_micro = usage_info.ru_utime.tv_usec + usage_info.ru_utime.tv_sec * 1000000;
    double user_time_sec = double(user_time_micro) / 1000000.0;
    logLine("    User time (seconds): " + std::to_string(user_time_sec));

    long long system_time_micro = usage_info.ru_stime.tv_usec + usage_info.ru_stime.tv_sec * 1000000;
    double system_time_sec = double(system_time_micro) / 1000000.0;
    logLine("    System time (seconds): " + std::to_string(system_time_sec));

    static std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
    long long wall_time_micro = std::chrono::duration_cast<std::chrono::microseconds>(now - _resourceLoggingStartPoint).count();
    double wall_time_sec = double(wall_time_micro) / 1000000.0;
    double percent_cpu = 100.0 * (user_time_sec + system_time_sec) / wall_time_sec;
    logLine("    Elapsed (wall clock) time (seconds): " + std::to_string(wall_time_sec));
    logLine("    Percent CPU this job got: " + std::to_string(percent_cpu) + "%");


    // missing: average (un)shared text size, stack size, total size --> all 0 for

    logLine("    Maximum resident set size (kbytes): " + std::to_string(usage_info.ru_maxrss));
    logLine("    Major (requiring I/O) page faults: " + std::to_string(usage_info.ru_majflt));
    logLine("    Minor (reclaiming a frame) page faults: " + std::to_string(usage_info.ru_minflt));

    logLine("    Voluntary context switches: " + std::to_string(usage_info.ru_nvcsw));
    logLine("    Involuntary context switches: " + std::to_string(usage_info.ru_nivcsw));

    logLine("    File system inputs: " + std::to_string(usage_info.ru_inblock));
    logLine("    File system outputs: " + std::to_string(usage_info.ru_oublock));

#else
    logLine("No resource usage information available. Only works under linux.")
#endif
}

void Logger::logCmdArgs(int argc, char **argv) {
    for(int i = 0; i < argc; i++) _cmdArgs.push_back(argv[i]);
}

void Logger::_logCmdArgs() {
    logLine("parameters:");
    std::string currLine = "    ";
    for (auto & arg : _cmdArgs) {
        if (arg.size() >= 2 && arg[0] == '-' && arg[1] == '-') {
            // argument is param name
            if (currLine.size() > 0) logLine(currLine);
            currLine = "    ";
        }

        currLine += arg + " ";
    }

    if (currLine.size() > 0) logLine(currLine);
}

void Logger::logError(std::exception_ptr eptr) {
    if (!eptr) {
        logLine("ERROR: std::bad_exception");
        throw std::bad_exception();
    } else {
        try { std::rethrow_exception(eptr); }
        catch (const std::exception &e) { logLine(std::string("ERROR: ") + e.what()); }
        catch (const std::string &e) { logLine(std::string("ERROR: ") + e); }
        catch (const char *e) { logLine(std::string("ERROR: ") + e); }
        catch (...) { logLine("ERROR: unknown error occurred."); }
    }
}
