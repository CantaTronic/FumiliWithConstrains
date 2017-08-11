
#pragma once

#include <vector>

class PDFGen {
protected:
  std::vector<float> x;
  virtual void PDF() = 0;
  virtual void Uniform() = 0;
  virtual void _Generate(const char * ofs_name, unsigned N, bool is_pdf);
public:
  PDFGen(unsigned nDim);
  virtual void Generate(const char * ofs_name, unsigned N) { _Generate(ofs_name, N, true); }
  virtual void GenerateUniform(const char * ofs_name, unsigned N) { _Generate(ofs_name, N, false); }
};
