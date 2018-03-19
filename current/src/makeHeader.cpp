/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "makeHeader.hpp"

Header::Header(std::string file_name, std::string H_gaurd) {
  header_file.open(file_name);
  H_gaurd = "AUTO_GEN_HEADER_" + H_gaurd;

  header_file << "// Auto Generated header file\n";
  header_file << "// Made from makeHeader.cpp\n\n";
  header_file << "#ifndef " << H_gaurd << "_H\n";
  header_file << "#define " << H_gaurd << "_H\n";
  header_file << "\n#include \"constants.hpp\"\n\n";
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
  end = true;
  // for each l in lines:
  for (auto &l : lines) l = "";
}

void Header::Set_RetrunType(std::string RetrunType) { r_type = RetrunType; }

void Header::Set_FuncName(std::string FuncName) { f_name = FuncName; }

void Header::Set_FuncInputs(std::string FuncInputs) { f_input = FuncInputs; }

void Header::Set_Function(TString Function) {
  func = std::string(Function.Data());
}

void Header::AddText(std::string TextAdd) { a_text = TextAdd; }

void Header::AddLine(std::string LineAdd) { lines.push_back(LineAdd); }
void Header::AddLine(std::string LineAdd, bool ending) {
  lines.push_back(LineAdd);
  end = ending;
}

void Header::AddComment(std::string CommAdd) { c_text = CommAdd; }

void Header::WriteFunction() {
  if (!r_type.empty() && !f_name.empty() && !f_input.empty() && !func.empty()) {
    header_file << r_type << " " << f_name << "(";
    header_file << f_input << ") {\n";

    if (c_text.length() > 0) header_file << "\t// " << c_text << "\n";
    if (a_text.length() > 0) header_file << "\t" << a_text << "\n";
    if (lines.size() > 0)
      for (auto &l : lines)
        if (!l.empty() && end) {
          header_file << "\t" << l << ";\n";
        } else {
          header_file << "\t" << l;
        }

    header_file << "\treturn " << func;
    header_file << ";\n}\n\n";
    Header::NewFunction();
  } else {
    std::cerr << "Cannot Write to header:" << std::endl
              << "r_type = " << !r_type.empty() << std::endl
              << "f_name = " << !f_name.empty() << std::endl
              << "f_input = " << !f_input.empty() << std::endl
              << "func = " << !func.empty() << std::endl;
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

  Header::AddLine("double a = " + std::to_string(a));
  Header::AddLine("double m = " + std::to_string(m));
  Header::AddLine("double s = " + std::to_string(s));
  Header::AddLine("double p = ((x - m) / s)");

  Header::Set_Function("(INV_SQRT_2PI / s * a * std::exp(-0.5f * p * p))");
  Header::WriteFunction();

  Header::NewFunction();
  Header::Set_RetrunType("bool");
  Header::Set_FuncName("between_" + name);
  Header::Set_FuncInputs("double x");
  Header::AddLine("double m = " + std::to_string(m));
  Header::AddLine("double s = " + std::to_string(s));
  Header::AddLine("bool above = (x <= m+s*N_SIGMA)");
  Header::AddLine("bool below = (x >= m-s*N_SIGMA)");
  Header::Set_Function("(above && below)");
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
