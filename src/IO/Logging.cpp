#include "Logging.h"

#include <utility>

Log::Logger::Logger(std::string prefix) : _prefix(std::move(prefix)), _nullstr() {
    _start_time = std::chrono::steady_clock::now();
    _nullstr.setstate(std::ios_base::badbit);
}

std::ostream &Log::Logger::error() {
    if (active_level >= Error ) {
        std::cout << _prefix;
        auto now = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else return _nullstr;
}

std::ostream &Log::Logger::info() {
    if (active_level >= Info) {
        std::cout << _prefix;
        auto now = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else return _nullstr; 
}

std::ostream &Log::Logger::warn() {
    if (active_level >= Warn) {
        std::cout << _prefix;
        auto now = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else return _nullstr; 
}

std::ostream &Log::Logger::debug() {
    if (active_level >= Debug) {
        std::cout << _prefix;
        auto now = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else return _nullstr; 
}

std::ostream &Log::Logger::trace() {
    if (active_level >= All) {
        std::cout << _prefix;
        auto now = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else return _nullstr; 
}
void Log::Logger::progress(const std::string &prefix, int step, int max) {
    if (step == 0) std::cout << _prefix << prefix;

    int digits = 1;
    int div = 10;
    while ((max / div) > 0) {
        digits++;
        div *= 10;
    }
    if (step > 0) for (int i = 0; i < digits; i++) std::cout << '\b';

    int step_digits = 1;
    div = 10;
    while ((step / div) > 0) {
        step_digits++;
        div *= 10;
    }

    for (int i = 0; i < digits - step_digits; i++) std::cout << ' ';
    std::cout << step;

    if (step >= max-1) {
        for (int i = 0; i < digits; i++)
            std::cout << '\b';
        for (int i = 0; i < prefix.size(); i++)
            std::cout << '\b';
        for (int i = 0; i < _prefix.size(); i++)
            std::cout << '\b';
    }
    //if (step >= max) std::cout << '\n';
    std::cout.flush();
}

void Log::init(Log::Level level) {
    active_level = level;
    general = std::make_unique<Logger>("[GENERAL] ");
    simulation = std::make_unique<Logger>("[SIMULATION] ");
    io = std::make_unique<Logger>("[IO] ");
}

void Log::finalize() {
    general = nullptr;
    simulation = nullptr;
    io = nullptr;
}

Log::Level Log::active_level = Info;
std::unique_ptr<Log::Logger> Log::general = nullptr;
std::unique_ptr<Log::Logger> Log::simulation = nullptr;
std::unique_ptr<Log::Logger> Log::io = nullptr;

const std::string Log::level_to_string(Log::Level level) {
    switch (level) {
    case Error:
        return "Error";
    case Info:
        return "Info";
    case Warn:
        return "Warn";
    case Debug:
        return "Debug";
    case All:
        return "All";
    }
    return "";
}
