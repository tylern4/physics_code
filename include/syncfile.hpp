#ifndef SYNCFILE_H_GUARD
#define SYNCFILE_H_GUARD

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include "csv_data.hpp"

// https://stackoverflow.com/questions/33596910/c-how-to-have-multiple-threads-write-to-a-file

class SyncFile {
 public:
  SyncFile(const std::string& path);
  ~SyncFile();
  bool write(const csv_data& data);
  bool write(const std::string& data);
  bool writeToFile();

 private:
  std::string _path;
  std::fstream _csv_output;
  std::mutex _writerMutex;
  std::queue<csv_data> _writeQueue;
};

#endif