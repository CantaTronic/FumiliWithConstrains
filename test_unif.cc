
#include "PDFGenSquare.h"
#include "MyPDF.h"

const int N = 100000; // кол-во чисел в файле
const float par[] = {0.5, 0.3, 0.8, 0.1}; // конст-ты pdf

int main() {
  PDFGenSquare pdf_g(MyPDF::nDim);
  pdf_g.SetParameters(MyPDF::nPar, par);
  float x[] = {1., 1.};
  pdf_g.SetPDF(MyPDF::PDF, MyPDF::PDF(x, par));
  pdf_g.Generate("unif_.dat", N);
  pdf_g.GenerateUniform("mc_.dat", N*5);
  return 0;
}
