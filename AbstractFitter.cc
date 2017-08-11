
#include "AbstractFitter.h"
#include <iostream>

void AbstractFumiliFitter::Minimize(const std::vector<FitParameter> & par) {
  // initialize fitter
  fumili = new TFumili(par.size());
  fumili->SetParNumber(par.size());
  // set parameters
  for(unsigned i = 0; i < par.size(); i++) {
    if(par[i].type == FitParameter::IS_FIXED) {
      fumili->SetParameter(i, par[i].name.c_str(),
                           par[i].val, par[i].step, 1., 1.);
    } else if(par[i].type == FitParameter::IS_FREE) {
      fumili->SetParameter(i, par[i].name.c_str(),
                           par[i].val, par[i].step, 0., 0.);
    } else {
      fumili->SetParameter(i, par[i].name.c_str(),
                           par[i].val, par[i].step,
                           par[i].min, par[i].max);
    }
  }
  // set objective function
  currentFitter = this;
  fumili->SetFCN(AbstractFitter::FCN);
  // init. sgz
  parameters = 0;
  objective_function = 0;
  gradients = 0;
  second_derivatives = fumili->GetZ();
  steps = fumili->GetPL0();
  // minimize
  double arglist[2] = {150, 0.01}; // TODO implement set n_iterations and tolerance
  fumili->ExecuteCommand("FUMILI", arglist, 2);
  // finalize
  fumili->PrintResults(1, 99);
  delete fumili;
  parameters = 0;
  objective_function = 0;
  gradients = 0;
  second_derivatives = 0;
  steps = 0;
}

// objective function
double AbstractFumiliLogLikelihoodFitter::FCN() {
  int n_par = fumili->GetNumberTotalParameters();
  int n_free_par = fumili->GetNumberFreeParameters();
  // fill SGZ with zeros
  objective_function = 0.;
  memset(gradients, 0, sizeof(double)*n_par);
  memset(second_derivatives, 0, sizeof(double)*n_free_par*(n_free_par+1)/2);
  // calculate SGZ
  for (unsigned iev = 0; iev < nev; iev++) {
    pdf = PDF(iev);
    objective_function += -log(pdf);
    // Derivatives
    double dpdf[n_par];
    for(int i = 0, iz = 0; i < n_par; i++) {
      if(steps[i] <= 0.) continue;
      // Calculation of gradient
      dpdf[i] = -dPDFdp(i, iev)/pdf;
      gradients[i] += dpdf[i];
      for(int j = 0; j <= i; j++) {
        if(steps[j] <= 0.) continue;
        second_derivatives[iz] += dpdf[i]*dpdf[j];
        iz++;
      }
    }
  }
  // add normalization
  norm = Norm();
  double dnorm[n_par];
  for(int i = 0; i < n_par; i++) {
    dnorm[i] = dNormdx(i)/norm;
  }
  objective_function += nev*log(norm);
  for(int i = 0, iz = 0; i < n_par; i++) {
    if(steps[i] <= 0.) continue;
    for(int j = 0; j <= i; j++) {
      if(steps[j] <= 0.) continue;
      second_derivatives[iz] +=
          dnorm[i]*gradients[j] + gradients[i]*dnorm[j] + nev*dnorm[i]*dnorm[j];
      iz++;
    }
  }
  for(int i = 0; i < n_par; i++) {
    if(steps[i] <= 0.) continue;
    gradients[i] += nev*dnorm[i];
  }
//   Print("S", 1, &objective_function);
//   Print("G", n_par, gradients);
//   Print("Z", (n_free_par+1)*n_free_par/2, second_derivatives);
  return objective_function;
}
double AbstractFumiliLogLikelihoodFitter::dPDFdp(unsigned i_par, unsigned i_ev) {
  double pval_save = parameters[i_par];
  parameters[i_par] += 1e-3;
  double ret = (PDF(i_ev)-pdf)*1e3;
  parameters[i_par] = pval_save;
  return ret;
}
double AbstractFumiliLogLikelihoodFitter::dNormdx(unsigned i_par) {
  double pval_save = parameters[i_par];
  parameters[i_par] += 1e-3;
  double ret = (Norm()-norm)*1e3;
  parameters[i_par] = pval_save;
  return ret;
}

AbstractFitter * AbstractFitter::currentFitter = 0;
// static objective function wrapper
void AbstractFitter::FCN(int & n_par, double * grad, double & ret_val, double * par, int flag) {
  currentFitter->SetParameters(par);
  currentFitter->SetGradients(grad);
  ret_val = currentFitter->FCN();
}
void AbstractFitter::Print(const char * name, int n, double * d) {
  std::cout<<name<<":";
  for(int i = 0; i < n; i++)
    std::cout<<" "<<d[i];
  std::cout<<std::endl;
}

void FitParameterSetter::AddParameter(std::string name, double val, double step) {
  parameters.resize(parameters.size()+1);
  parameters.back().name = name;
  parameters.back().val = val;
  parameters.back().step = step;
  parameters.back().type = FitParameter::IS_FREE;
}
void FitParameterSetter::AddParameter(std::string name, double val, double min, double max, double step) {
  parameters.resize(parameters.size()+1);
  parameters.back().name = name;
  parameters.back().val = val;
  parameters.back().step = step;
  parameters.back().type = FitParameter::IS_LIMITED;
  parameters.back().min = min;
  parameters.back().max = max;
}
void FitParameterSetter::AddFixedParameter(std::string name, double val, double step) {
  parameters.resize(parameters.size()+1);
  parameters.back().name = name;
  parameters.back().val = val;
  parameters.back().step = step;
  parameters.back().type = FitParameter::IS_FIXED;
}
