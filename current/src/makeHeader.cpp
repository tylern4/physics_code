/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "makeHeader.hpp"

Header::Header(std::string file_name) {
  header_file.open(file_name);
  std::string name_name = "AUTO_GEN_HEADER";

  header_file << "// Auto Generated header file\n";
  header_file << "// Made from makeHeader.cpp\n\n";
  header_file << "#ifndef " << name_name << "_H\n";
  header_file << "#define " << name_name << "_H\n";
  header_file << "\n";
}

Header::~Header() {
  header_file << "#endif\n" << std::endl;
  header_file.close();
}

void Header::NewFunction() {
  r_type = "";
  f_name = "";
  f_input = "";
  func = "";
  a_text = "";
  c_text = "";
  // for each l in lines:
  for (auto &l : lines)
    l = "";
}

void Header::Set_RetrunType(std::string RetrunType) { r_type = RetrunType; }

void Header::Set_FuncName(std::string FuncName) { f_name = FuncName; }

void Header::Set_FuncInputs(std::string FuncInputs) { f_input = FuncInputs; }

void Header::Set_Function(TString Function) {
  func = std::string(Function.Data());
}

void Header::AddText(std::string TextAdd) { a_text = TextAdd; }

void Header::AddLine(std::string LineAdd) { lines.push_back(LineAdd); }

void Header::AddComment(std::string CommAdd) { c_text = CommAdd; }

void Header::WriteFunction() {
  if (!r_type.empty() && !f_name.empty() && !f_input.empty() && !func.empty()) {
    header_file << r_type << " " << f_name << "(";
    header_file << f_input << ") {\n";

    if (c_text.length() > 0)
      header_file << "\t// " << c_text << "\n";
    if (a_text.length() > 0)
      header_file << "\t" << a_text << "\n";
    if (lines.size() > 0)
      for (auto &l : lines)
        header_file << "\t" << l << ";\n";

    header_file << "\treturn " << func;
    header_file << ";\n}\n\n";
    Header::NewFunction();
  } else {
    std::cerr << "Cannot Write to header:" << std::endl;
  }
}

void Header::WriteGaussian(std::string name, double a, double m, double s) {
  /*
  double normal_pdf(double x, double m, double s){
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
  }
  */
  Header::NewFunction();
  Header::Set_RetrunType("double");
  Header::Set_FuncName("gaussian_" + name);
  Header::Set_FuncInputs("double x");

  Header::AddLine("static const float inv_sqrt_2pi = 0.3989422804014327");
  Header::AddLine("a = " + std::to_string(a));
  Header::AddLine("m = " + std::to_string(m));
  Header::AddLine("s = " + std::to_string(s));
  Header::AddLine("double p = (x - m) / s");

  Header::Set_Function("inv_sqrt_2pi / s * a * std::exp(-0.5f * p * p)");
  Header::WriteFunction();
}

void Header::MakeFunction(std::string RetrunType, std::string FuncName,
                          std::string FuncInputs, std::string Function) {
  Header::NewFunction();
  Header::Set_RetrunType(RetrunType);
  Header::Set_FuncName(FuncName);
  Header::Set_FuncInputs(FuncInputs);
  Header::Set_Function(Function);
  Header::WriteFunction();
}
