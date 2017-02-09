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
  header_file << "\n\n#endif\n" << std::endl;
  header_file.close();
}

void Header::Set_RetrunType(std::string RetrunType) { r_type = RetrunType; }

void Header::Set_FuncName(std::string FuncName) { f_name = FuncName; }

void Header::Set_FuncInputs(std::string FuncInputs) { f_input = FuncInputs; }

void Header::Set_Function(std::string Function) { func = Function; }

void Header::AddText(std::string TextAdd) { a_text = TextAdd; }

void Header::AddComment(std::string CommAdd) { c_text = CommAdd; }

void Header::WriteFunction() {
  std::cout << "header_file" << std::endl;
  header_file << r_type << " " << f_name << "(";
  header_file << f_input << ") {\n";
  if (c_text.length() > 0)
    header_file << "\t// " << c_text << "\n";
  if (a_text.length() > 0)
    header_file << "\t" << a_text << "\n";
  header_file << "\treturn " << func;
  header_file << ";\n}";
}

void Header::MakeFunction(std::string RetrunType, std::string FuncName,
                          std::string FuncInputs, std::string Function) {
  Header::Set_RetrunType(RetrunType);
  Header::Set_FuncName(FuncName);
  Header::Set_FuncInputs(FuncInputs);
  Header::Set_Function(Function);
  Header::WriteFunction();
}
