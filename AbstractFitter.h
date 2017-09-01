
#pragma once

#include <vector>
#include <string>

class FitParameters {
protected:
  unsigned n_par, n_free_par;
  std::vector<std::string> names;
  std::vector<double> params;
  enum Type {IS_FIXED, IS_FREE, IS_LIMITED};
  std::vector<Type> types;
  std::vector<double> steps;
  std::vector<double> minima;
  std::vector<double> maxima;
public:
  FitParameters(): n_par(0), n_free_par(0) { }
  virtual void AddParameter(std::string name, double val, double step = 0.01);
  virtual void AddParameter(std::string name, double val, double min, double max, double step = 0.01);
  virtual void AddFixedParameter(std::string name, double val);
  virtual void SetParameters(const double * par);
  virtual unsigned NPar() { return n_par; }
  virtual unsigned NFreePar() { return n_free_par; }
  virtual const char * Name(unsigned i) { return names[i].c_str(); }
  virtual double Val(unsigned i) { return params[i]; }
  virtual double Step(unsigned i) { return steps[i]; }
  virtual double Min(unsigned i) { return minima[i]; }
  virtual double Max(unsigned i) { return maxima[i]; }
  virtual bool IsFixed(unsigned i) { return types[i] == IS_FIXED; }
};

class ObjectiveFunction: public FitParameters {
protected:
  double last_val;
public:
  unsigned fcn_call;
  virtual void Eval(double * par, double & val, double * grad = 0, double * Z = 0) = 0;
  virtual double LastVal() { return last_val; }
  virtual ~ObjectiveFunction() { }
  static void Print(const char * name, unsigned n, double * d);
  void PrintParameters();
};

class LogLikelihoodObjFunc: public ObjectiveFunction {
protected:
  unsigned nev;
  double pdf, norm;
//   unsigned i_ev, i_par;
  virtual double PDF(unsigned i_ev) = 0;
  virtual double dPDFdp(unsigned i_ev, unsigned i_par);
  virtual double Norm() { return 1.; }
  virtual double dNormdx(unsigned i_par);
public:
  virtual void Eval(double * _par, double & val, double * grad = 0, double * Z = 0);
  virtual ~LogLikelihoodObjFunc() { }
};

class Minimizer {
protected:
  ObjectiveFunction * objFunc;
public:
  Minimizer(): objFunc(0) { }
  Minimizer(ObjectiveFunction * f): objFunc(f) { }
  virtual ~Minimizer() { }
  virtual void SetObjFunc(ObjectiveFunction * f) { objFunc = f; }
  virtual void Minimize() = 0;
};

class TVirtualFitter;

class ROOTMinimizer: public Minimizer {
protected:
  static void FCN(int & n_par, double * gradients, double & ret_val, double * par, int flag);
  static TVirtualFitter * gFitter;
  static ObjectiveFunction * gObjFunc;
public:
  enum Strategy {FUMILI, SIMPLEX, MIGRAD, MIGRAD_G};
  static Strategy gStrategy;
  ROOTMinimizer() { }
  ROOTMinimizer(Strategy s) { gStrategy = s; }
  ROOTMinimizer(ObjectiveFunction * f): Minimizer(f) { }
  ROOTMinimizer(Strategy s, ObjectiveFunction * f): Minimizer(f) { gStrategy = s; }
  virtual ~ROOTMinimizer() { }
  virtual void Minimize();
};
