
#include "Minimizer.h"
#include "ObjectiveFunction.h"
#include <TFumili.h>
#include <TFitter.h>
#include <TMinuit.h>

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
