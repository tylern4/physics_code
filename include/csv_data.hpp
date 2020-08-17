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
  float e_p;
  float e_cx;
  float e_cy;
  float e_cz;
  float pip_p;
  float pip_cx;
  float pip_cy;
  float pip_cz;
  int helicty;
  std::string type;
  size_t hash;

  friend std::ostream &operator<<(std::ostream &os, const csv_data &d) {
    os << d.electron_sector << ",";
    os << d.w << ",";
    os << d.q2 << ",";
    os << d.theta << ",";
    os << d.phi << ",";
    os << d.mm2 << ",";
    os << d.e_p << ",";
    os << d.e_cx << ",";
    os << d.e_cy << ",";
    os << d.e_cz << ",";
    os << d.pip_p << ",";
    os << d.pip_cx << ",";
    os << d.pip_cy << ",";
    os << d.pip_cz << ",";
    os << d.helicty << ",";
    os << d.type << ",";
    os << d.hash;
    return os;
  }
};

#endif