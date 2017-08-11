
#include "TestFitter.h"
#include "MyPDF.h"
#include <iostream>
#include <fstream>

int main() {
  FitParameterSetter par;
  par.AddFixedParameter("#alpha1", .5);
  par.AddParameter("#beta1", 0);
  par.AddParameter("#alpha2", 1);
  par.AddParameter("#beta2", 0);
  TestFumiliFitter fitter;
  fitter.ReadData("unif_.dat");
  fitter.ReadNormData("mc_.dat");
  fitter.Minimize(par);
}

bool TestFumiliFitter::ReadData(const std::string filename) {
  nev = 0;
  float data_i[MyPDF::nDim];
  std::ifstream ifs(filename.c_str());
  if(!ifs || ifs.eof()) return false;
  while(true) {
    for(unsigned iDim = 0; iDim < MyPDF::nDim; iDim++)
      ifs >> data_i[iDim];
    if(!ifs || ifs.eof()) break;
    nev++;
    for(unsigned iDim = 0; iDim < MyPDF::nDim; iDim++)
      data.push_back(data_i[iDim]);
  }
  ifs.close();
  std::cout<<"TestFumiliFitter::ReadData("<<filename<<"): read "<<nev<<" events"<<std::endl;
  return nev;
}
bool TestFumiliFitter::ReadNormData(const std::string filename) {
  nev_norm = 0;
  float data_i[MyPDF::nDim];
  std::ifstream ifs(filename.c_str());
  if(!ifs || ifs.eof()) return false;
  while(true) {
    for(unsigned iDim = 0; iDim < MyPDF::nDim; iDim++)
      ifs >> data_i[iDim];
    if(!ifs || ifs.eof()) break;
    nev_norm++;
    for(unsigned iDim = 0; iDim < MyPDF::nDim; iDim++)
      data_norm.push_back(data_i[iDim]);
  }
  ifs.close();
  std::cout<<"TestFumiliFitter::ReadNormData("<<filename<<"): read "<<nev_norm<<" events"<<std::endl;
  return nev_norm;
}
double TestFumiliFitter::PDF(unsigned i_ev) {
  return MyPDF::PDF(&data.data()[i_ev*MyPDF::nDim], parameters);
}
double TestFumiliFitter::Norm() {
  double norm = 0;
  for(unsigned iev = 0; iev < nev_norm; iev++) {
    norm += MyPDF::PDF(&data_norm.data()[iev*MyPDF::nDim], parameters);
  }
  norm /= nev_norm;
  return norm;
}
