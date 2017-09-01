
#pragma once

class ObjectiveFunction;

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
