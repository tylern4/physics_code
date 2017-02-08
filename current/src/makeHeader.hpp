/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include <iostream>
#include <fstream>
#include <string>

#ifndef MAKE_HEADER_H
#define MAKE_HEADER_H

class Header {
private:
  std::ofstream header_file;

public:
  Header(std::string file_name);
  ~Header();
};

#endif
