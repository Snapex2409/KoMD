#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

namespace Log {

enum Level {
    Error,
    Info,
    Warn,
    Debug,
    All
};

[[maybe_unused]] const std::string level_to_string(Level level);

extern Level active_level;

class Logger {
  public:
    explicit Logger(std::string prefix);

    /// error > info > warn > debug > trace
    std::ostream& error();
    /// error > info > warn > debug > trace
    std::ostream& info();
    /// error > info > warn > debug > trace
    std::ostream& warn();
    /// error > info > warn > debug > trace
    std::ostream& debug();
    /// error > info > warn > debug > trace
    std::ostream& trace();
    /// progress bar
    void progress(const std::string& prefix, int step, int max);

  private:
    /// starting time
    std::chrono::steady_clock::time_point _start_time;
    /// what the logger should always print first
    const std::string _prefix;
    /// empty stream to write nothing
    std::ofstream _nullstr;
};

/**
 * Creates all loggers and sets the level
 * */
void init(Level level = Info);

/**
 * Frees all memory
 * */
void finalize();

/// logger with prefix [GENERAL]
extern std::unique_ptr<Logger> general;
/// logger with prefix [SIMULATION]
extern std::unique_ptr<Logger> simulation;
/// logger with prefix [IO]
extern std::unique_ptr<Logger> io;

}; // namespace Log
