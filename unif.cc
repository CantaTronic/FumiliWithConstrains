
#include "unif.h"
#include <TRandom.h>
#include <iomanip>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>

PDFGen::PDFGen(const char * ofs_name, const char * mc_name)
    : pdf(PDFGen::pdf_default),
      x_min(0), x_max(1), y_min(0), y_max(1), pdf_max(0),
      rootfile(0), hist_xy(0), hist_x(0), hist_y(0) {
  gRandom->SetSeed(0);
  ofs.open(ofs_name);
  if(mc_name) mc.open(mc_name);
}
void PDFGen::SetLimits(float _x_min, float _x_max, float _y_min, float _y_max) {
  x_min = _x_min;
  x_max = _x_max;
  y_min = _y_min;
  y_max = _y_max;
}
void PDFGen::SetPDF(float (*_pdf)(float x, float y), float _pdf_max) {
  pdf = _pdf;
  pdf_max = _pdf_max;
}
void PDFGen::InitROOT(const char * rootfile_name) {
  rootfile = new TFile(rootfile_name, "RECREATE");
  hist_xy = new TH2F("hist_xy", "Bar-chart for x & y", 10, x_min, x_max, 10, x_min, x_max);
  hist_x = new TH1F("hist_x", "Bar-chart for x", 10, x_min, x_max);
  hist_y = new TH1F("hist_y", "Bar-chart for y", 10, x_min, x_max);
  hist_x->SetMinimum(0);
  hist_y->SetMinimum(0);
}
void PDFGen::Generate(unsigned N) {
  unsigned i = 0;
  while (i < N) {
    float x = gRandom->Uniform(x_min, x_max);
    float y = gRandom->Uniform(x_min, x_max);
    float z = gRandom->Uniform(0, pdf_max);
    if (z < pdf(x, y)) {
      ofs<<std::setw(15)<<x<<std::setw(15)<<y<<std::endl;
      if (rootfile) {
        hist_x->Fill(x);
        hist_y->Fill(y);
        hist_xy->Fill(x,y);
      }
      i++;
    }
  }
  if(mc.is_open()) {
    for(unsigned i = 0; i < N*10; i++) {
      float x = gRandom->Uniform(x_min, x_max);
      float y = gRandom->Uniform(x_min, x_max);
      mc<<std::setw(15)<<x<<std::setw(15)<<y<<std::endl;
    }
  }
}
void PDFGen::DrawROOT(TF1 * fx, TF1 * fy, TF2 * fxy) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  TCanvas * c = new TCanvas("my_canv", "The canvas", 1000, 700);
  if (fx) hist_x->Fit(fx, "R");
  hist_x->Draw();
  c->SaveAs("X_hist.pdf");
  if (fy) hist_y->Fit(fy, "R");
  hist_y->Draw();
  c->SaveAs("Y_hist.pdf");
  if (fxy) hist_xy->Fit(fxy, "R");
  hist_xy->Draw("col");
  c->SaveAs("XY_hist.pdf");
  c->Close();
}
PDFGen::~PDFGen() {
  if (rootfile) {
    rootfile->Write();
    rootfile->Close();
  }
  ofs.close();
  if(mc.is_open()) mc.close();
}
float PDFGen::pdf_default(float x, float y) {
  return 1;
}
