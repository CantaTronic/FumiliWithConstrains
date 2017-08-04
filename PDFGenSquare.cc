
#include "PDFGenSquare.h"
#include <TRandom.h>
#include <iostream>

PDFGenSquare::PDFGenSquare(unsigned nDim)
: PDFGen(nDim), pdf(pdf_default) {
  min.assign(nDim, 0);
  max.assign(nDim, 1);
}
void PDFGenSquare::SetLimits(unsigned dim, float _min, float _max) {
  if(dim >= x.size()) {
    std::cerr<<"PDFGenSquare::SetLimits:"<<std::endl
             <<"  requested dimension "<<dim<<" is more than number of dimensions"<<x.size()<<std::endl;
    return;
  }
  min[dim] = _min;
  max[dim] = _max;
}
void PDFGenSquare::SetParameters(unsigned n, const float * par) {
  p.assign(par, par+n);
}
void PDFGenSquare::SetPDF(float (*_pdf)(const float *, const float *), float _pdf_max) {
  pdf = _pdf;
  pdf_max = _pdf_max;
}
void PDFGenSquare::PDF() {
  while(true) {
    for(unsigned i = 0; i < x.size(); i++)
      x[i] = gRandom->Uniform(min[i], max[i]);
    float z = gRandom->Uniform(0, pdf_max);
    if(z < pdf(x.data(), p.data())) return;
  }
}
void PDFGenSquare::Uniform() {
  for(unsigned i = 0; i < x.size(); i++)
    x[i] = gRandom->Uniform(min[i], max[i]);
}
float PDFGenSquare::pdf_default(const float * x, const float * p) {
  return 1;
}
