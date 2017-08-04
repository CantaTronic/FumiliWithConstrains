
#pragma once

#include "PDFGen.h"

class PDFGenSquare: public PDFGen {
private:
  std::vector<float> p;
  std::vector<float> min;
  std::vector<float> max;
  float (*pdf)(const float * x, const float * p);
  float pdf_max;
  static float pdf_default(const float * x, const float * p);
  virtual void PDF();
  virtual void Uniform();
public:
  PDFGenSquare(unsigned nDim);
  virtual void SetLimits(unsigned dim, float _min, float _max);
  virtual void SetParameters(unsigned n, const float * p);
  virtual void SetPDF(float (*_pdf)(const float * x, const float * p), float _pdf_max);
};
