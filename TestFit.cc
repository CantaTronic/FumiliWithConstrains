
#include "TestFit.h"
#include "MyPDF.h"
#include <iostream>
#include <fstream>

int main() {
  TestFit fcn;
  fcn.ReadData("unif_.dat");
  fcn.ReadNormData("mc_.dat");
  fcn.Minimize(ROOTMinimizer::MIGRAD, "MIGRAD GRAD_NONE");
  fcn.Minimize(ROOTMinimizer::MIGRAD_G, "MIGRAD GRAD_FORCE");
  fcn.Minimize(ROOTMinimizer::SIMPLEX, "SIMPLEX");
  fcn.Minimize(ROOTMinimizer::FUMILI, "FUMILI");
}

const double TestFit::defaultParams[] = {.5, 0., 1., 0.};
TestFit::TestFit() {
  AddFixedParameter("#alpha1", defaultParams[0]);
  AddParameter("#beta1", defaultParams[1]);
  AddParameter("#alpha2", defaultParams[2]);
  AddParameter("#beta2", defaultParams[3]);
}
void TestFit::Minimize(ROOTMinimizer::Strategy strategy, const char * title) {
  SetParameters(defaultParams);
  std::cout<<"\n ==== "<<title<<" ===="<<std::endl;
  timer.Start();
  ROOTMinimizer * minimizer = new ROOTMinimizer(strategy, this);
  minimizer->Minimize();
  delete minimizer;
  timer.Stop();
  timer.Print();
  timer.Reset();
}
bool TestFit::ReadData(const std::string filename) {
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
  std::cout<<"TestFit::ReadData("<<filename<<"): read "<<nev<<" events"<<std::endl;
  return nev;
}
bool TestFit::ReadNormData(const std::string filename) {
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
  std::cout<<"TestFit::ReadNormData("<<filename<<"): read "<<nev_norm<<" events"<<std::endl;
  return nev_norm;
}
double TestFit::PDF(unsigned i_ev) {
  return MyPDF::PDF(&data.data()[i_ev*MyPDF::nDim], params.data());
}
double TestFit::Norm() {
  double norm = 0;
  for(unsigned iev = 0; iev < nev_norm; iev++) {
    norm += MyPDF::PDF(&data_norm.data()[iev*MyPDF::nDim], params.data());
  }
  norm /= nev_norm;
  return norm;
}
