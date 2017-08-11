
#include "PDFGen.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <TRandom.h>

PDFGen::PDFGen(unsigned nDim) {
  x.resize(nDim);
  gRandom->SetSeed(0);
}

void PDFGen::_Generate(const char * ofs_name, unsigned N, bool is_pdf) {
  std::ofstream ofs(ofs_name);
  if(!ofs || !ofs.is_open()) {
    std::cerr<<"Cannot open file "<<ofs_name<<" for writing!"<<std::endl;
    return;
  }
  for (unsigned i = 0; i < N; i++) {
    if(is_pdf) PDF();
    else Uniform();
    for(unsigned j = 0; j < x.size(); j++)
      ofs<<std::setw(15)<<x[j];
    ofs<<std::endl;
  }
  ofs.close();
}
