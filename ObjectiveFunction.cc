
#include "ObjectiveFunction.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>

void LogLikelihoodObjFunc::Eval(double * par, double & val, double * grad, double * Z) {
  fcn_call++;
  SetParameters(par);
//   Print("P", n_par, par);
  // fill SGZ with zeros
  val = 0.;
  if(grad) memset(grad, 0, sizeof(double)*n_par);
  if(Z) memset(Z, 0, sizeof(double)*n_free_par*(n_free_par+1)/2);
  // calculate SGZ
  for(unsigned iev = 0; iev < nev; iev++) {
    pdf = PDF(iev);
    val += -log(pdf);
    if(!grad) continue;
    // Derivatives
    double dpdf[n_par];
    for(unsigned i_par = 0, i_z = 0; i_par < n_par; i_par++) {
      if(IsFixed(i_par)) continue;
      // Calculation of gradient
      dpdf[i_par] = -dPDFdp(iev, i_par)/pdf;
      grad[i_par] += dpdf[i_par];
      if(!Z) continue;
      for(unsigned j_par = 0; j_par <= i_par; j_par++) {
        if(IsFixed(j_par)) continue;
        Z[i_z] += dpdf[i_par]*dpdf[j_par];
        i_z++;
      }
    }
  }
  // add normalization
  norm = Norm();
  val += nev*log(norm);
  last_val = val;
//   Print("S", 1, &val);
  if(!grad) return;
  double dnorm[n_par];
  for(unsigned i_par = 0; i_par < n_par; i_par++) {
    if(IsFixed(i_par)) continue;
    dnorm[i_par] = dNormdx(i_par)/norm;
  }
  if(Z) {
    for(unsigned i_par = 0, i_z = 0; i_par < n_par; i_par++) {
      if(IsFixed(i_par)) continue;
      for(unsigned j_par = 0; j_par <= i_par; j_par++) {
        if(IsFixed(j_par)) continue;
        Z[i_z] += dnorm[i_par]*grad[j_par]
                + grad[i_par]*dnorm[j_par]
                + nev*dnorm[i_par]*dnorm[j_par];
        i_z++;
      }
    }
//   Print("Z", (n_free_par+1)*n_free_par/2, Z);
  }
  for(unsigned i_par = 0; i_par < n_par; i_par++) {
    if(IsFixed(i_par)) continue;
    grad[i_par] += nev*dnorm[i_par];
  }
//   Print("G", n_par, grad);
//   PrintParameters();
}
double LogLikelihoodObjFunc::dPDFdp(unsigned i_ev, unsigned i_par) {
  double pval_save = params[i_par];
  params[i_par] += 1e-3;
  double ret = (PDF(i_ev)-pdf)*1e3;
  params[i_par] = pval_save;
  return ret;
}
double LogLikelihoodObjFunc::dNormdx(unsigned i_par) {
  double pval_save = params[i_par];
  params[i_par] += 1e-3;
  double ret = (Norm()-norm)*1e3;
  params[i_par] = pval_save;
  return ret;
}

void ObjectiveFunction::Print(const char * name, unsigned n, double * d) {
  std::cout<<name<<":";
  for(unsigned i = 0; i < n; i++)
    std::cout<<" "<<d[i];
  std::cout<<std::endl;
}
void ObjectiveFunction::PrintParameters() {
  std::cout<<"FCN call"<<std::setw(5)<<fcn_call<<", parameters";
  for(unsigned i = 0; i < NPar(); i++)
    std::cout<<std::setw(13)<<params[i];
/*
  if(gStrategy == MIGRAD || gStrategy == SIMPLEX) {
    std::cout<<' '<<static_cast<TFitter*>(root_fitter)->GetMinuit()->fCfrom;
    std::cout<<' '<<static_cast<TFitter*>(root_fitter)->GetMinuit()->fCstatu;
  }
*/
  std::cout<<std::endl;
}

void FitParameters::AddParameter(std::string name, double val, double step) {
  n_par++;
  n_free_par++;
  names.push_back(name);
  params.push_back(val);
  steps.push_back(step);
  types.push_back(IS_FREE);
  minima.push_back(0.);
  maxima.push_back(0.);
}
void FitParameters::AddParameter(std::string name, double val, double min, double max, double step) {
  n_par++;
  n_free_par++;
  names.push_back(name);
  params.push_back(val);
  steps.push_back(step);
  types.push_back(IS_LIMITED);
  minima.push_back(min);
  maxima.push_back(max);
}
void FitParameters::AddFixedParameter(std::string name, double val) {
  n_par++;
  names.push_back(name);
  params.push_back(val);
  steps.push_back(0.01);
  types.push_back(IS_FIXED);
  minima.push_back(0.);
  maxima.push_back(0.);
}
void FitParameters::SetParameters(const double * par) {
  for(unsigned i = 0; i < params.size(); i++) {
    if(types[i] == IS_FIXED) continue;
    params[i] = par[i];
  }
}
