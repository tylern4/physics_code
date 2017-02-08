/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "makeHeader.hpp"

Header::Header(std::string file_name) {
  header_file.open(file_name);
  std::string name_name = "HEADER";

  header_file << "// Auto Generated header file\n";
  header_file << "// Made from makeHeader.cpp\n\n";
  header_file << "#ifndef " << name_name << "_H\n";
  header_file << "#define " << name_name << "_H\n";
  header_file << "\n\n";
}

Header::~Header() {
  header_file << "#endif\n" << std::endl;
  header_file.close();
}

Header::PutFunction(std::string FuncName, std::string FuncInputs,
                    std::string Function) {}
