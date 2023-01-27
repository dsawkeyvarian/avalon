#include <iostream>
#include <fmt/core.h>
#include "logger.h"

#define SIZEOFARRAY(arr) (sizeof(arr)/sizeof(*arr))

static const char* log_level_names[] = {
  "DEBUG",
  "INFO",
  "WARNING",
  "ERROR"
};

static_assert(static_cast<size_t>(LogLevel::NUM_LOG_LEVELS) == SIZEOFARRAY(log_level_names));

ILogger& ILogger::logWithLevel(LogLevel level, const std::string &s)
{
    const std::string msg = fmt::format("[{}]: {}\n", 
      log_level_names[static_cast<int>(level)], s);
    this->log(msg);
    return *this;
}
ILogger& ILogger::debug(const std::string &s) {
    return logWithLevel(LogLevel::DEBUG, s);
}
ILogger& ILogger::info(const std::string &s) {
    return logWithLevel(LogLevel::INFO, s);
}
ILogger& ILogger::warning(const std::string &s) {
    return logWithLevel(LogLevel::WARNING, s);
}
ILogger& ILogger::error(const std::string &s) {
    return logWithLevel(LogLevel::ERROR, s);
}

COutAndFileLogger::COutAndFileLogger() {}
COutAndFileLogger::COutAndFileLogger(const fs::path& file) {
    m_of.open(file);
    if(!m_of.good())
        throw std::runtime_error(fmt::format("Unable to open file {} for logging", file.string()));
}
COutAndFileLogger::~COutAndFileLogger() {
    if(m_of.is_open())
        m_of.close();
}
ILogger& COutAndFileLogger::log(const std::string &msg) {
    std::cout << msg;
    if(m_of.is_open())
        m_of << msg;
    return *this;
}

class Logger
{
public:
  Logger(const std::string& filename, bool to_cout) 
    : m_to_cout(to_cout) 
  {
    m_of.open(filename);
  }

  Logger()
    : m_to_cout(true)
  {}

  ~Logger() {
    //if(m_of.is_open())
    //  m_of.close();
  }

  Logger(const Logger &) = delete;
  Logger& logWithLevel(LogLevel level, const std::string &s)
  {
    const std::string msg = fmt::format("[{}]: {}\n", 
      log_level_names[static_cast<int>(level)], s);
    return log(msg);
  }
  Logger& log(const std::string &msg) {
      if(m_of.is_open())
      m_of << msg;
    if(m_to_cout)
      std::cout << msg;
    return *this;
  }
  Logger& debug(const std::string &s) {
    return logWithLevel(LogLevel::DEBUG, s);
  }
  Logger& info(const std::string &s) {
    return logWithLevel(LogLevel::INFO, s);
  }
  Logger& warning(const std::string &s) {
    return logWithLevel(LogLevel::WARNING, s);
  }
  Logger& error(const std::string &s) {
    return logWithLevel(LogLevel::ERROR, s);
  }

private:
  bool m_to_cout;
  std::ofstream m_of;
};