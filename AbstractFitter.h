
#pragma once

#include <vector>
#include <string>
#include <TFumili.h>
#include <TFitter.h>

struct FitParameter {
  std::string name;
  double val;
  enum {IS_FIXED, IS_FREE, IS_LIMITED} type;
  double step;
  double min;
  double max;
};

class FitParameterSetter {
protected:
  std::vector<FitParameter> parameters;
public:
  virtual void AddParameter(std::string name, double val, double step = 0.01);
  virtual void AddParameter(std::string name, double val, double min, double max, double step = 0.01);
  virtual void AddFixedParameter(std::string name, double val);
  virtual operator std::vector<FitParameter>&() { return parameters; }
};

class AbstractFitter {
protected:
  unsigned fcn_call;
  TVirtualFitter * root_fitter;
  double * parameters;
  std::vector<bool> is_par_fixed;
  double objective_function;
  double * gradients;
  double * second_derivatives;
  static AbstractFitter * currentFitter;
  static void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
  virtual void InitFitter(const std::vector<FitParameter> & par);
  virtual void MinimizeFumili(const std::vector<FitParameter> & par);
  virtual void MinimizeMinuit(const std::vector<FitParameter> & par);
  virtual double FCNFumili() = 0;
  virtual double FCNMigrad() = 0;
public:
  enum Strategy {FUMILI, MIGRAD, SIMPLEX} strategy;
  enum {GRAD_NONE, GRAD_CHECK, GRAD_FORCE} do_user_gradients;
  AbstractFitter(Strategy _strategy = MIGRAD)
  : fcn_call(0), root_fitter(0), parameters(0),
    objective_function(0), gradients(0), second_derivatives(0),
    strategy(_strategy), do_user_gradients(GRAD_CHECK) { }
  virtual ~AbstractFitter() { }
  virtual bool HasGradients() { return gradients; }
  virtual void SetGradients(double * grad) { gradients = grad; }
  virtual void SetParameters(double * par) { parameters = par; }
  virtual void Minimize(const std::vector<FitParameter> & par);
  virtual double FCN();
  static void Print(const char * name, int n, double * d);
  virtual void PrintParameters();
};

class AbstractLogLikelihoodFitter: public AbstractFitter {
protected:
  unsigned nev;
  double pdf, norm;
  virtual double PDF(unsigned i_ev) = 0;
  virtual double dPDFdp(unsigned i_par, unsigned i_ev);
  virtual double Norm() { return 1.; }
  virtual double dNormdx(unsigned i_par);
  virtual double FCNFumili();
  virtual double FCNMigrad();
public:
  virtual ~AbstractLogLikelihoodFitter() { }
};
