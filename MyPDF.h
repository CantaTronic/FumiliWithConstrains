
#pragma once

namespace MyPDF {
  extern const unsigned nDim;
  extern const unsigned nPar;
  inline double PDF(const double * x, const double * p) {
    return 1 + p[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[1] + p[3]*x[1]*x[1];
  }
  inline float PDF(const float * x, const float * p) {
    return 1 + p[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[1] + p[3]*x[1]*x[1];
  }
}