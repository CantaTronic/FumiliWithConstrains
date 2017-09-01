
#pragma once

#include "ObjectiveFunction.h"
#include "Minimizer.h"
#include <TStopwatch.h>

class TestFit: public LogLikelihoodObjFunc {
protected:
  static const double defaultParams[];
  TStopwatch timer;
  std::vector<double> data, data_norm;
  unsigned nev_norm;
public:
  TestFit();
  virtual ~TestFit() { }
  virtual void Minimize(ROOTMinimizer::Strategy strategy, const char * title);
  virtual bool ReadData(const std::string filename);
  virtual bool ReadNormData(const std::string filename);
  virtual double PDF(unsigned i_ev);
  virtual double Norm();
};
