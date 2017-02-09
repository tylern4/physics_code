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
  std::string r_type = "";
  std::string f_name = "";
  std::string f_input = "";
  std::string func = "";
  std::string a_text = "";
  std::string c_text = "";

public:
  Header(std::string file_name);
  ~Header();
  void Set_RetrunType(std::string RetrunType);
  void Set_FuncName(std::string FuncName);
  void Set_FuncInputs(std::string FuncInputs);
  void Set_Function(std::string Function);
  void AddText(std::string TextAdd);
  void AddComment(std::string CommAdd);
  void MakeFunction(std::string RetrunType, std::string FuncName,
                    std::string FuncInputs, std::string Function);
  void WriteFunction();
};

#endif
