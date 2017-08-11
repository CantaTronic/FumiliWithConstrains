
#pragma once

#include <vector>
#include <string>
#include <TFumili.h>

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
  virtual void AddFixedParameter(std::string name, double val, double step = 0.01);
  virtual operator std::vector<FitParameter>&() { return parameters; }
};

class AbstractFitter {
protected:
  double * parameters;
  double objective_function;
  double * gradients;
  static AbstractFitter * currentFitter;
  static void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
public:
  AbstractFitter(): parameters(0), objective_function(0), gradients(0) { }
  virtual bool HasGradients() { return gradients; }
  virtual void SetGradients(double * grad) { gradients = grad; }
  virtual void SetParameters(double * par) { parameters = par; }
  virtual void Minimize(const std::vector<FitParameter> & par) = 0;
  virtual double FCN() = 0;
  static void Print(const char * name, int n, double * d);
};

class AbstractFumiliFitter: public AbstractFitter {
protected:
  TFumili * fumili;
  double * second_derivatives;
  double * steps;
public:
  AbstractFumiliFitter(): fumili(0), second_derivatives(0), steps(0) { }
  virtual void Minimize(const std::vector<FitParameter> & par);
  virtual double FCN() = 0;
};

class AbstractFumiliLogLikelihoodFitter: public AbstractFumiliFitter {
protected:
  unsigned nev;
  double pdf, norm;
  virtual double PDF(unsigned i_ev) = 0;
  virtual double dPDFdp(unsigned i_par, unsigned i_ev);
  virtual double Norm() { return 1.; }
  virtual double dNormdx(unsigned i_par);
public:
  virtual double FCN();
};
