#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  int electron_sector;
  float w;
  float q2;
  float theta;
  float phi;
  float mm2;
  int helicty;
  std::string type;
  size_t hash;

  friend std::ostream &operator<<(std::ostream &os, const csv_data &d) {
    os << std::setprecision(10);
    os << d.electron_sector << ",";
    os << d.w << ",";
    os << d.q2 << ",";
    os << d.theta << ",";
    os << d.phi << ",";
    os << d.mm2 << ",";
    os << std::to_string(d.helicty) << ",";
    os << d.type;
    // << ",";
    // os << std::hex << d.hash;
    return os;
  }
};

#endif