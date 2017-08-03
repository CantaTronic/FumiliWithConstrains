
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
}

int main(int arc, char ** argv) {
  FCN::ReadData("unif.dat");
  // initialize fitter
  FCN::fumili = new TFumili;
  FCN::fumili->SetParNumber(2);
  // set parameters
  FCN::fumili->SetParameter(0, // number
                       "#alpha1", // name
                       0.5, // start value
                       0.01, // start step
                       0.5, // lower limit
                       0.5); // upper limit
  FCN::fumili->SetParameter(1, "#alpha2", 1, 0.01, 0, 0);
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
  Print("P", 2, par);
  Print("PL0", 2, PL0);
//   exit(0);
  ret_val = 0.;
  for(int i = 0; i < n_par; i++) {
    gradients[i] = .0;
    for(int j = 0; j <= i; j++)
      Z[(i+1)*i/2+j] = .0;
  }
  for (int iev = 0; iev < nev; iev++) {
    double pdf_unnormalized = 1. + par[0]*data[iev][0] + par[1]*data[iev][1];
    double normalization   = 1. + 0.5*par[0] + 0.5*par[1];
    double pdf = pdf_unnormalized/normalization;
    ret_val -= log(pdf);
    // Derivatives
    double df[2];
    for(int i = 0; i < n_par; i++) {
      double dch = data[iev][i];
      double dzn = 0.5;
      double ys = (dch - pdf*dzn)/normalization;
      df[i] = ys/pdf;
    }
    int iz = 0;
    for(int i = 0; i < n_par; i++) {
      if(PL0[i] <= .0) continue;
      // Calculation of gradient
      gradients[i] = gradients[i] - df[i];
      for(int j = 0; j <= i; j++) {
        if(PL0[j] <= .0) continue;
        Z[iz] = Z[iz] + df[i]*df[j];
        iz++;
      }
    }
  }
  Print("S", 1, &ret_val);
  Print("G", 2, gradients);
  Print("Z", 3, Z);
}
// read data from file
void FCN::ReadData(const char * filename) {
  std::ifstream ifs("unif.dat");
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
