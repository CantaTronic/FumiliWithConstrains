
#include "AbstractFitter.h"
#include <TMinuit.h>
#include <iostream>
#include <iomanip>

void AbstractFitter::Minimize(const std::vector<FitParameter> & par) {
  // init.
  fcn_call = 0;
  is_par_fixed.clear();
  for(unsigned i = 0; i < par.size(); i++)
    is_par_fixed.push_back(par[i].type == FitParameter::IS_FIXED);
  parameters = 0;
  objective_function = 0;
  gradients = 0;
  second_derivatives = 0;
  switch(strategy) {
    case FUMILI:
      MinimizeFumili(par); break;
    case MIGRAD: case SIMPLEX:
      MinimizeMinuit(par); break;
  }
  // finalize
  delete root_fitter;
  fcn_call = 0;
  parameters = 0;
  objective_function = 0;
  gradients = 0;
  second_derivatives = 0;
  is_par_fixed.clear();
}
void AbstractFitter::MinimizeFumili(const std::vector<FitParameter> & par) {
  // initialize fitter
  TFumili * fumili = new TFumili(par.size());
  fumili->SetParNumber(par.size());
  root_fitter = fumili;
  InitFitter(par);
  // init. z
  second_derivatives = fumili->GetZ();
  // minimize
  double arglist[2] = {150, 0.01}; // TODO implement set n_iterations and tolerance
  fumili->ExecuteCommand("FUMILI", arglist, 2);
  // finalize
  fumili->PrintResults(1, objective_function);
}
void AbstractFitter::MinimizeMinuit(const std::vector<FitParameter> & par) {
  // initialize fitter
  TFitter * minuit = new TFitter(par.size());
  double arglist[2];
  arglist[0] = -1;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  root_fitter = minuit;
  InitFitter(par);
  if(strategy == MIGRAD) {
    arglist[0] = 1;
    if(do_user_gradients == GRAD_NONE)
      minuit->ExecuteCommand("SET NOGRADIENT", arglist, 0);
    else if(do_user_gradients == GRAD_FORCE)
      minuit->ExecuteCommand("SET GRADIENT", arglist, 1);
    else
      minuit->ExecuteCommand("SET GRADIENT", arglist, 0);
  }
//   arglist[0] = 0;
//   root_fitter->ExecuteCommand("SET STRATEGY", arglist, 1);
  // minimize
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  // Minimize
  // TODO implement set n_iterations and tolerance
  arglist[0] = 5000; // number of function calls
  arglist[1] = 0.01; // tolerance
  if(strategy == MIGRAD)
    minuit->ExecuteCommand("MIGRAD", arglist, 2);
  if(strategy == SIMPLEX)
    minuit->ExecuteCommand("SIMPLEX", arglist, 2);
}
double AbstractFitter::FCN() {
  fcn_call++;
  switch(strategy) {
    case FUMILI:
      return FCNFumili();
    case MIGRAD: case SIMPLEX:
      return FCNMigrad();
  }
  return 0;
}

// objective function
double AbstractLogLikelihoodFitter::FCNFumili() {
  int n_par = root_fitter->GetNumberTotalParameters();
  int n_free_par = root_fitter->GetNumberFreeParameters();
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
      if(is_par_fixed[i]) continue;
      // Calculation of gradient
      dpdf[i] = -dPDFdp(i, iev)/pdf;
      gradients[i] += dpdf[i];
      for(int j = 0; j <= i; j++) {
        if(is_par_fixed[j]) continue;
        second_derivatives[iz] += dpdf[i]*dpdf[j];
        iz++;
      }
    }
  }
  // add normalization
  norm = Norm();
  double dnorm[n_par];
  for(int i = 0; i < n_par; i++) {
    if(is_par_fixed[i]) continue;
    dnorm[i] = dNormdx(i)/norm;
  }
  objective_function += nev*log(norm);
  for(int i = 0, iz = 0; i < n_par; i++) {
    if(is_par_fixed[i]) continue;
    for(int j = 0; j <= i; j++) {
      if(is_par_fixed[j]) continue;
      second_derivatives[iz] +=
          dnorm[i]*gradients[j] + gradients[i]*dnorm[j] + nev*dnorm[i]*dnorm[j];
      iz++;
    }
  }
  for(int i = 0; i < n_par; i++) {
    if(is_par_fixed[i]) continue;
    gradients[i] += nev*dnorm[i];
  }
//   Print("P", n_par, parameters);
//   Print("S", 1, &objective_function);
//   Print("G", n_par, gradients);
//   Print("Z", (n_free_par+1)*n_free_par/2, second_derivatives);
//   PrintParameters();
  return objective_function;
}
// objective function
double AbstractLogLikelihoodFitter::FCNMigrad() {
  bool do_gradients = do_user_gradients;
  if(strategy == SIMPLEX)
    do_gradients = false;
  if(static_cast<TFitter*>(root_fitter)->GetMinuit()->fCfrom == "HESSE   ")
    do_gradients = false;
  int n_par = root_fitter->GetNumberTotalParameters();
  // fill SG with zeros
  objective_function = 0.;
  if(do_gradients)
    memset(gradients, 0, sizeof(double)*n_par);
  // calculate SG
  for (unsigned iev = 0; iev < nev; iev++) {
    pdf = PDF(iev);
    objective_function += -log(pdf);
    if(!do_gradients) continue;
    // Derivatives
    for(int i = 0; i < n_par; i++) {
      if(is_par_fixed[i]) continue;
      // Calculation of gradient
      gradients[i] += -dPDFdp(i, iev)/pdf;
    }
  }
  // add normalization
  norm = Norm();
  objective_function += nev*log(norm);
  if(do_gradients) {
    for(int i = 0; i < n_par; i++) {
      if(is_par_fixed[i]) continue;
      gradients[i] += nev*dNormdx(i)/norm;
    }
  }
//   PrintParameters();
  return objective_function;
}
double AbstractLogLikelihoodFitter::dPDFdp(unsigned i_par, unsigned i_ev) {
  double pval_save = parameters[i_par];
  parameters[i_par] += 1e-3;
  double ret = (PDF(i_ev)-pdf)*1e3;
  parameters[i_par] = pval_save;
  return ret;
}
double AbstractLogLikelihoodFitter::dNormdx(unsigned i_par) {
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
void AbstractFitter::InitFitter(const std::vector<FitParameter> & par) {
  // set parameters
  for(unsigned i = 0; i < par.size(); i++) {
    root_fitter->SetParameter(i, par[i].name.c_str(),
                              par[i].val, par[i].step,
                              par[i].min, par[i].max);
    if(par[i].type == FitParameter::IS_FIXED) root_fitter->FixParameter(i);
  }
  // set objective function
  currentFitter = this;
  root_fitter->SetFCN(FCN);
}
void AbstractFitter::Print(const char * name, int n, double * d) {
  std::cout<<name<<":";
  for(int i = 0; i < n; i++)
    std::cout<<" "<<d[i];
  std::cout<<std::endl;
}
void AbstractFitter::PrintParameters() {
  std::cout<<"FCN call"<<std::setw(5)<<fcn_call<<", parameters";
  for(int i = 0; i < root_fitter->GetNumberTotalParameters(); i++)
    std::cout<<std::setw(13)<<parameters[i];
  if(strategy == MIGRAD || strategy == SIMPLEX) {
    std::cout<<' '<<static_cast<TFitter*>(root_fitter)->GetMinuit()->fCfrom;
    std::cout<<' '<<static_cast<TFitter*>(root_fitter)->GetMinuit()->fCstatu;
  }
  std::cout<<std::endl;
}

void FitParameterSetter::AddParameter(std::string name, double val, double step) {
  parameters.resize(parameters.size()+1);
  parameters.back().name = name;
  parameters.back().val = val;
  parameters.back().step = step;
  parameters.back().type = FitParameter::IS_FREE;
  parameters.back().min = 0.;
  parameters.back().max = 0.;
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
void FitParameterSetter::AddFixedParameter(std::string name, double val) {
  parameters.resize(parameters.size()+1);
  parameters.back().name = name;
  parameters.back().val = val;
  parameters.back().step = 0.01;
  parameters.back().type = FitParameter::IS_FIXED;
  parameters.back().min = 0.;
  parameters.back().max = 0.;
}
