
#include <TFumili.h>
#include <fstream>
#include <iostream>
#include <cmath>

namespace FCN {
  const int nev = 100000;
  double data[nev][2];
  void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
  void ReadData(const char * filename);
  TFumili * fumili;
  void Print(const char * name, int n, double * d);
  double PDF(double * x, double * p);
  double Norm(double * p);
  double dPDFdx(int np, double * x);
  double dNormdx(int np);
}

int main(int arc, char ** argv) {
  FCN::ReadData("unif_.dat");
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
  double * Z = fumili->GetZ();
  double * PL0 = fumili->GetPL0();
  std::cout<<"N: "<<n_par<<std::endl;
  Print("P", n_par, par);
  Print("PL0", n_par, PL0);
//   exit(0);
  ret_val = 0.;
  int iz = 0;
  for(int i = 0; i < n_par; i++) {
    gradients[i] = .0;
    if(PL0[i] <= .0) continue;
    for(int j = 0; j <= i; j++) {
      if(PL0[j] <= .0) continue;
      Z[iz] = .0;
      iz++;
    }
  }
  for (int iev = 0; iev < nev; iev++) {
//     double pdf_unnormalized = 1. + par[0]*data[iev][0] + par[1]*data[iev][1];
//     double normalization   = 1. + 0.5*par[0] + 0.5*par[1];
    double pdf = PDF(data[iev], par);
    double norm = Norm(par);
    ret_val += -log(pdf/norm);
    // Derivatives
    double df[n_par];
    for(int i = 0; i < n_par; i++) {
      df[i] = -dPDFdx(i, data[iev])/pdf + dNormdx(i)/norm;
    }
    int iz = 0;
    for(int i = 0; i < n_par; i++) {
      if(PL0[i] <= .0) continue;
      // Calculation of gradient
      gradients[i] += df[i];
      for(int j = 0; j <= i; j++) {
        if(PL0[j] <= .0) continue;
        Z[iz] += df[i]*df[j];
        iz++;
      }
    }
  }
  Print("S", 1, &ret_val);
  Print("G", n_par, gradients);
  Print("Z", (n_par+1)*n_par/2, Z);
}
// read data from file
void FCN::ReadData(const char * filename) {
  std::ifstream ifs(filename);
  for(int i = 0; i < nev; i++)
    ifs >> data[i][0] >> data[i][1];
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
  return 1 + p[0]/2. + p[1]/3. + p[2]/2. + p[3]/3.;
}
double FCN::dPDFdx(int np, double * x) {
  double p[4] = {0.};
  p[np] = 1;
  return p[0]*x[0] + p[1]*x[0]*x[0] + p[2]*x[1] + p[3]*x[1]*x[1];
}
double FCN::dNormdx(int np) {
  double p[4] = {0.};
  p[np] = 1;
  return p[0]/2. + p[1]/3. + p[2]/2. + p[3]/3.;
}
