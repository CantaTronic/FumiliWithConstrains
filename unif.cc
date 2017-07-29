
#include <TRandom.h>
#include <fstream>
#include <iomanip>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;

const int N = 100000; //кол-во чисел в файле
const float x_min = 0, x_max = 1; //границы отрезка
const float alpha1 = 0.5, alpha2 = 0.8, pdf_max = 1 + (alpha1 + alpha2)*x_max; //конст-ты pdf
const float pdf1_max = 1 + (alpha1)*x_max;

float pdf(float x, float y) {
  return 1 + alpha1*x + alpha2*y;
}

float pdf1(float x) {
  return 1 + alpha1*x;
}

int main() {
  gRandom->SetSeed(0);
  ofstream ofs("Unif.dat"); //файл для записи чисел
  TFile * rootfile = new TFile("unif.root", "recreate");
  TH2F * hist_xy = new TH2F("hist_xy", "Bar-chart for x & y", 10, x_min, x_max, 10, x_min, x_max);
  TH1F * hist_x = new TH1F("hist_x", "Bar-chart for x", 10, x_min, x_max);
  TH1F * hist_y = new TH1F("hist_y", "Bar-chart for y", 10, x_min, x_max);
  hist_x->SetMinimum(0);
  hist_y->SetMinimum(0);
  int i=0;
  while (i < N) {
    float x = gRandom->Uniform(x_min,x_max);
    float y = gRandom->Uniform(x_min,x_max);
    float z = gRandom->Uniform(0,pdf_max);
    if (pdf(x,y) > z) {
//       ofs<<setw(15)<<x<<setw(15)<<y<<endl;
      hist_x->Fill(x);
      hist_y->Fill(y);
      hist_xy->Fill(x,y);
      i++;
    }
  }
  ofs.close();

  TF1 * f1 = new TF1("f1", "[0]*(1+[1]*x)", 0, 1);
  TF2 * f2 = new TF2("f2", "[0]*(1+[1]*x+[2]*y)", 0, 1, 0, 1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  TCanvas * c = new TCanvas("my_canv", "The canvas", 1000, 700);
  hist_x->Fit(f1, "R");
  hist_x->Draw();
  c->SaveAs("X_hist.pdf");
  hist_y->Fit(f1, "R");
  hist_y->Draw();
  c->SaveAs("Y_hist.pdf");
  hist_xy->Fit(f2, "R");
  hist_xy->Draw("col");
  c->SaveAs("XY_hist.pdf");
  c->Close();

  rootfile->Write();
  rootfile->Close();

  return 0;
}
