
#include "unif.h"
#include <TF1.h>
#include <TF2.h>

const int N = 1000000; // кол-во чисел в файле
const float alpha1 = 0.5, beta1 = 0.3, alpha2 = 0.8, beta2 = 0.1; // конст-ты pdf
// const float pdf1_max = 1 + (alpha1)*x_max;

float f_pdf(float x, float y) {
  return 1 + alpha1*x + beta1*x*x + alpha2*y + beta2*y*y;
}

// float pdf1(float x) {
//   return 1 + alpha1*x;
// }

int main() {
  PDFGen pdf_g("unif_.dat"/*, "mc_.dat"*/);
  pdf_g.SetPDF(f_pdf, 1+alpha1+beta1+alpha2+beta2/*f_pdf(1, 1)*/);
  pdf_g.InitROOT("unif_.root");
  pdf_g.Generate(N);
//   TF1 * fx = new TF1("fx", "[0]*(1+[1]*x+[2]*x*x)", 0, 1);
//   TF1 * fy = new TF1("fy", "[0]*(1+[1]*x)", 0, 1);
//   TF2 * f2 = new TF2("fxy", "[0]*(1+[1]*x+[2]*x*x+[3]*y+[4]*y*y)", 0, 1, 0, 1);
  pdf_g.DrawROOT(0, 0/*, f2*/);
  return 0;
}
