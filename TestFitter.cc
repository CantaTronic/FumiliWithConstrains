
#include "TestFitter.h"
#include "MyPDF.h"
#include <TStopwatch.h>
#include <iostream>
#include <fstream>

int main() {
  FitParameterSetter par;
  par.AddFixedParameter("#alpha1", .5);
  par.AddParameter("#beta1", 0);
  par.AddParameter("#alpha2", 1);
  par.AddParameter("#beta2", 0);
  TestFitter * fitter = new TestFitter;
  fitter->ReadData("unif_.dat");
  fitter->ReadNormData("mc_.dat");
  TStopwatch timer;
  std::cout<<"\n ==== MIGRAD GRAD_NONE ===="<<std::endl;
  timer.Start();
  fitter->strategy = AbstractFitter::MIGRAD;
  fitter->do_user_gradients = AbstractFitter::GRAD_NONE;
  fitter->Minimize(par);
  timer.Stop();
  timer.Print();
  timer.Reset();
  std::cout<<"\n ==== MIGRAD GRAD_CHECK ===="<<std::endl;
  timer.Start();
  fitter->strategy = AbstractFitter::MIGRAD;
  fitter->do_user_gradients = AbstractFitter::GRAD_CHECK;
  fitter->Minimize(par);
  timer.Stop();
  timer.Print();
  timer.Reset();
  std::cout<<"\n ==== MIGRAD GRAD_FORCE ===="<<std::endl;
  timer.Start();
  fitter->strategy = AbstractFitter::MIGRAD;
  fitter->do_user_gradients = AbstractFitter::GRAD_FORCE;
  fitter->Minimize(par);
  timer.Stop();
  timer.Print();
  timer.Reset();
  std::cout<<"\n ==== SIMPLEX ===="<<std::endl;
  timer.Start();
  fitter->strategy = AbstractFitter::SIMPLEX;
  fitter->Minimize(par);
  timer.Stop();
  timer.Print();
  timer.Reset();
  std::cout<<"\n ==== FUMILI ===="<<std::endl;
  timer.Start();
  fitter->strategy = AbstractFitter::FUMILI;
  fitter->Minimize(par);
  timer.Stop();
  timer.Print();
  timer.Reset();
  delete fitter;
}

bool TestFitter::ReadData(const std::string filename) {
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
  std::cout<<"TestFitter::ReadData("<<filename<<"): read "<<nev<<" events"<<std::endl;
  return nev;
}
bool TestFitter::ReadNormData(const std::string filename) {
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
  std::cout<<"TestFitter::ReadNormData("<<filename<<"): read "<<nev_norm<<" events"<<std::endl;
  return nev_norm;
}
double TestFitter::PDF(unsigned i_ev) {
  return MyPDF::PDF(&data.data()[i_ev*MyPDF::nDim], parameters);
}
double TestFitter::Norm() {
  double norm = 0;
  for(unsigned iev = 0; iev < nev_norm; iev++) {
    norm += MyPDF::PDF(&data_norm.data()[iev*MyPDF::nDim], parameters);
  }
  norm /= nev_norm;
  return norm;
}
