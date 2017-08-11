
#pragma once

#include "AbstractFitter.h"

class TestFumiliFitter: public AbstractFumiliLogLikelihoodFitter {
  protected:
    std::vector<double> data, data_norm;
    unsigned nev_norm;
  public:
    virtual bool ReadData(const std::string filename);
    virtual bool ReadNormData(const std::string filename);
    virtual double PDF(unsigned i_ev);
    virtual double Norm();
};
