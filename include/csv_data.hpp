#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
 public:
  int electron_sector;
  float w;
  float q2;
  float theta;
  float phi;
  float mm2;
  int helicty;
  std::string type;

  static std::string header() { return "electron_sector,w,q2,theta,phi,mm2,helicty,type"; }

  friend std::ostream &operator<<(std::ostream &os, const csv_data &d) {
    os << std::setprecision(10);
    os << d.electron_sector << ",";
    os << d.w << ",";
    os << d.q2 << ",";
    os << d.theta << ",";
    os << d.phi << ",";
    os << d.mm2 << ",";
    os << std::to_string(d.helicty) << ",";
    os << d.type << "\n";
    return os;
  }
};

struct fid_data : csv_data {
 public:
  int sector;
  float e_p;
  float e_theta;
  float e_phi;
  float theta;
  float phi;
  float x;
  float y;
  std::string type;

  static std::string header() { return "sector,e_p,e_theta,e_phi,theta,phi,x,y,type\n"; }

  const std::string print() {
    std::string out = "";
    out += std::to_string(sector) + ",";
    out += std::to_string(e_p) + ",";
    out += std::to_string(e_theta) + ",";
    out += std::to_string(e_phi) + ",";
    out += std::to_string(theta) + ",";
    out += std::to_string(phi) + ",";
    out += std::to_string(x) + ",";
    out += std::to_string(y) + ",";
    out += type + "\n ";
    return out;
  }
};

#endif