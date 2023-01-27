#pragma once
#include <filesystem>
namespace fs = std::filesystem;

fs::path getVirtualLinacRoot();
fs::path getOutputFolder();