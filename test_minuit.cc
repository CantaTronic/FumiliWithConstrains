
#include <TVirtualFitter.h>
#include <fstream>
#include <cmath>

namespace FCN {
  const int nev = 100000;
  double data[nev][2];
  void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
  void ReadData(const char * filename);
}

int main(int arc, char ** argv) {
  FCN::ReadData("unif.dat");
  // initialize fitter
  double arglist[2];
  delete TVirtualFitter::GetFitter();
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * minuit = TVirtualFitter::Fitter(0 /* zero */, 2 /* num. of param~s */);
  arglist[0] = -1; // no printing
  minuit->ExecuteCommand("SET PRINT", // command
                         arglist, // arguments
                         1); // number of arg~s
  // set parameters
  minuit->SetParameter(0, // number
                       "#alpha1", // name
                       1, // start value
                       0.01, // start step
                       -10, // lower limit
                       10); // upper limit
  minuit->SetParameter(1, "#alpha2", 1, 0.01, -10, 10);
  // set objective function
  minuit->SetFCN(FCN::FCN);
  // basic printing
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  // minimize
  arglist[0] = 5000; // number of function calls
  arglist[1] = 0.01; // tolerance
  minuit->ExecuteCommand("MIGRAD"/* minimizer used */, arglist, 2);
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
  ret_val = 0.;
  for (int iev = 0; iev < nev; iev++) {
    double pdf_unnormalized = 1. + par[0]*data[iev][0] + par[1]*data[iev][1];
    double normalization   = 1. + 0.5*par[0] + 0.5*par[1];
    double pdf = pdf_unnormalized/normalization;
    ret_val -= log(pdf);
  }
}
// read data from file
void FCN::ReadData(const char * filename) {
  std::ifstream ifs("unif.dat");
  for(int i = 0; i < nev; i++)
    ifs >> data[i][0] >> data[i][1];
  ifs.close();
}
