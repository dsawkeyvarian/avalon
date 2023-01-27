#pragma once
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

enum class LogLevel {
  DEBUG,
  INFO,
  WARNING,
  ERROR,
  NUM_LOG_LEVELS
};

class ILogger
{
public:
  ILogger() {}
  virtual ~ILogger() {}

  ILogger(const ILogger &) = delete;
  virtual ILogger& log(const std::string &msg) = 0;

  ILogger& logWithLevel(LogLevel level, const std::string &s);
  ILogger& debug(const std::string &s);
  ILogger& info(const std::string &s);
  ILogger& warning(const std::string &s);
  ILogger& error(const std::string &s);
};

class COutAndFileLogger : public ILogger {
public:
  COutAndFileLogger();
  COutAndFileLogger(const fs::path& file);
  ~COutAndFileLogger();

  ILogger& log(const std::string &msg) override;
private:
  std::ofstream m_of;
};

class SinkLogger : public ILogger {
public:
  ILogger& log(const std::string&) override {
    return *this;
  }
};