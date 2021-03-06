
#include "MyPDF.h"
#include <TFumili.h>
#include <fstream>
#include <iostream>
#include <cmath>

namespace FCN {
  const int nev = 100000;
  double data[nev][2];
  double mc[nev*5][2];
  void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
  void ReadData(const char * filename);
  void ReadNormData(const char * filename);
  TFumili * fumili;
  void Print(const char * name, int n, double * d);
  double PDF(double * x, double * p);
  double Norm(double * p);
  double dPDFdx(int np, double * x, double * p = 0, double pdf = 0);
  double dNormdx(int np, double * p = 0, double norm = 0);
}

int main(int arc, char ** argv) {
//   std::cout<<"Norm MC: "<<norm<<", ";
//   norm = FCN::Norm(p);
//   std::cout<<"Norm analytical: "<<norm<<std::endl;
  // start analysis
  FCN::ReadData("unif_.dat");
  FCN::ReadNormData("mc_.dat");
  // initialize fitter
  FCN::fumili = new TFumili;
  FCN::fumili->SetParNumber(4);
  // set parameters
  FCN::fumili->SetParameter(0, "#alpha1", .5, 0.01, .5, .5);
  FCN::fumili->SetParameter(1, "#beta1", 0, 0.01, 0, 0);
  FCN::fumili->SetParameter(2, "#alpha2", 1, 0.01, 0, 0);
  FCN::fumili->SetParameter(3, "#beta2", 0, 0.01, 0, 0);
  // set objective function
  FCN::fumili->SetFCN(FCN::FCN);
  // minimize
  FCN::fumili->Minimize();
  FCN::fumili->PrintResults(1, 99);
/*
  // get parameters
  double par[2], err[2];
  for(int i = 0; i < minuit->GetNumberTotalParameters(); i++) {
    if(par) par[i] = minuit->GetParameter(i);
    if(err) err[i] = minuit->GetParError(i);
  }
*/
  return 0;
}
// objective function
void FCN::FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag) {
  n_par = fumili->GetNumberTotalParameters();
  int n_free_par = fumili->GetNumberFreeParameters();
  double * Z = fumili->GetZ();
  double * PL0 = fumili->GetPL0();
//   std::cout<<"N: "<<n_par<<std::endl;
//   Print("P", n_par, par);
//   Print("PL0", n_par, PL0);
  // fill SGZ with zeros
  ret_val = 0.;
  memset(gradients, 0, sizeof(double)*n_par);
  memset(Z, 0, sizeof(double)*n_free_par*(n_free_par+1)/2);
  // calculate SGZ
  for (int iev = 0; iev < nev; iev++) {
    double pdf = PDF(data[iev], par);
    ret_val += -log(pdf);
    // Derivatives
    double dpdf[n_par];
    for(int i = 0, iz = 0; i < n_par; i++) {
      if(PL0[i] <= .0) continue;
      // Calculation of gradient
      dpdf[i] = -dPDFdx(i, data[iev], par, pdf)/pdf;
      gradients[i] += dpdf[i];
      for(int j = 0; j <= i; j++) {
        if(PL0[j] <= .0) continue;
        Z[iz] += dpdf[i]*dpdf[j];
        iz++;
      }
    }
  }
  // add normalization
  double norm = Norm(par);
  double dnorm[n_par];
  for(int i = 0; i < n_par; i++) {
    dnorm[i] = dNormdx(i, par, norm)/norm;
  }
  ret_val += nev*log(norm);
  for(int i = 0, iz = 0; i < n_par; i++) {
    if(PL0[i] <= .0) continue;
    for(int j = 0; j <= i; j++) {
      if(PL0[j] <= .0) continue;
      Z[iz] += dnorm[i]*gradients[j] + gradients[i]*dnorm[j] + nev*dnorm[i]*dnorm[j];
      iz++;
    }
  }
  for(int i = 0; i < n_par; i++) {
    if(PL0[i] <= .0) continue;
    gradients[i] += nev*dnorm[i];
  }
//   Print("S", 1, &ret_val);
//   Print("G", n_par, gradients);
//   Print("Z", (n_free_par+1)*n_free_par/2, Z);
}
// read data from file
void FCN::ReadData(const char * filename) {
  std::ifstream ifs(filename);
  for(int i = 0; i < nev; i++)
    ifs >> data[i][0] >> data[i][1];
  ifs.close();
}
void FCN::ReadNormData(const char * filename) {
  // check MC normalization
  std::ifstream ifs(filename);
  for(int i = 0; i < nev*5; i++)
    ifs >> mc[i][0] >> mc[i][1];
  ifs.close();
}
void FCN::Print(const char * name, int n, double * d) {
  std::cout<<name<<":";
  for(int i = 0; i < n; i++)
    std::cout<<" "<<d[i];
  std::cout<<std::endl;
}
double FCN::PDF(double * x, double * p) {
  return 1 + p[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[1] + p[3]*x[1]*x[1];
}
double FCN::Norm(double * p) {
//   return 1 + p[0]/2. + p[1]/3. + p[2]/2. + p[3]/3.;
  double norm = 0;
  for(int iev = 0; iev < nev*5; iev++) {
    norm += PDF(mc[iev], p);
  }
  norm /= nev*5;
  return norm;
}
double FCN::dPDFdx(int np, double * x, double * p, double pdf) {
  double _p[4] = {0.};
  if(p) { // numeric calculation
    memcpy(_p, p, 4*sizeof(double));
    _p[np] += 1e-3;
    double _pdf = PDF(x, _p);
    return (_pdf-pdf)*1e3;
  }
  _p[np] = 1;
  return _p[0]*x[0] + _p[1]*x[0]*x[0] + _p[2]*x[1] + _p[3]*x[1]*x[1];
}
double FCN::dNormdx(int np, double * p, double norm) {
  double _p[4] = {0.};
  if(p) { // numeric calculation
    memcpy(_p, p, 4*sizeof(double));
    _p[np] += 1e-3;
    double _norm = Norm(_p);
    return (_norm-norm)*1e3;
  }
  _p[np] = 1;
  return _p[0]/2. + _p[1]/3. + _p[2]/2. + _p[3]/3.;
}
