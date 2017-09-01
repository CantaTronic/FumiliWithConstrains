
#include "AbstractFitter.h"
#include <TFumili.h>
#include <TFitter.h>
#include <TMinuit.h>
#include <iostream>
#include <iomanip>

TVirtualFitter * ROOTMinimizer::gFitter = 0;
ObjectiveFunction * ROOTMinimizer::gObjFunc = 0;
ROOTMinimizer::Strategy ROOTMinimizer::gStrategy = ROOTMinimizer::MIGRAD;
// static objective function wrapper
void ROOTMinimizer::FCN(int & n_par, double * grad, double & ret_val, double * par, int flag) {
  switch(gStrategy) {
    case FUMILI:
      gObjFunc->Eval(par, ret_val, grad, static_cast<TFumili*>(gFitter)->GetZ());
      break;
    case SIMPLEX: case MIGRAD:
      gObjFunc->Eval(par, ret_val);
      break;
    case MIGRAD_G:
      if(static_cast<TFitter*>(gFitter)->GetMinuit()->fCfrom == "HESSE   ")
        gObjFunc->Eval(par, ret_val);
      else
        gObjFunc->Eval(par, ret_val, grad);
      break;
  }
}
void ROOTMinimizer::Minimize() {
  if(!objFunc) return;
  objFunc->fcn_call = 0;
  gObjFunc = objFunc;
  // initialize fitter
  double arglist[2];
  if(gStrategy == FUMILI) {
    gFitter = new TFumili(objFunc->NPar());
    static_cast<TFumili*>(gFitter)->SetParNumber(objFunc->NPar());
  }
  if(gStrategy == SIMPLEX || gStrategy == MIGRAD || gStrategy == MIGRAD_G) {
    gFitter = new TFitter(objFunc->NPar());
    arglist[0] = -1;
    gFitter->ExecuteCommand("SET PRINT", arglist, 1);
  }
  // set parameters
  for(unsigned i = 0; i < objFunc->NPar(); i++) {
    gFitter->SetParameter(i, objFunc->Name(i),
                          objFunc->Val(i), objFunc->Step(i),
                          objFunc->Min(i), objFunc->Max(i));
    if(objFunc->IsFixed(i)) gFitter->FixParameter(i);
  }
  // set objective function
  gFitter->SetFCN(FCN);
  // minimize
  if(gStrategy == FUMILI) {
    // TODO implement set n_iterations and tolerance
    arglist[0] = 150; // number of function calls
    arglist[1] = 0.01; // tolerance
    gFitter->ExecuteCommand("FUMILI", arglist, 2);
    gFitter->PrintResults(1, objFunc->LastVal());
  }
  if(gStrategy == SIMPLEX || gStrategy == MIGRAD || gStrategy == MIGRAD_G) {
    if(gStrategy == MIGRAD) {
      gFitter->ExecuteCommand("SET NOGRADIENT", arglist, 0);
    }
    if(gStrategy == MIGRAD_G) {
      arglist[0] = 1;
      gFitter->ExecuteCommand("SET GRADIENT", arglist, 1);
    }
//     arglist[0] = 0;
//     root_fitter->ExecuteCommand("SET STRATEGY", arglist, 1);
    arglist[0] = 0;
    gFitter->ExecuteCommand("SET PRINT", arglist, 1);
    // TODO implement set n_iterations and tolerance
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    if(gStrategy == SIMPLEX) {
      gFitter->ExecuteCommand("SIMPLEX", arglist, 2);
    }
    if(gStrategy == MIGRAD || gStrategy == MIGRAD_G) {
      gFitter->ExecuteCommand("MIGRAD", arglist, 2);
    }
  }
  // finalize
  delete gFitter;
}

// objective function
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
