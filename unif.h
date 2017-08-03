
#pragma once

#include <fstream>

class TFile;
class TH2F;
class TH1F;
class TF1;
class TF2;

class PDFGen {
private:
  std::ofstream ofs;
  float (*pdf)(float x, float y);
  float x_min, x_max, y_min, y_max;
  float pdf_max;
  TFile * rootfile;
  TH2F * hist_xy;
  TH1F * hist_x;
  TH1F * hist_y;
  static float pdf_default(float x, float y);
public:
  PDFGen(const char * ofs_name);
  void SetLimits(float _x_min, float _x_max, float _y_min, float _y_max);
  void SetPDF(float (*_pdf)(float x, float y), float _pdf_max);
  void InitROOT(const char * rootfile_name);
  void Generate(unsigned N);
  void DrawROOT(TF1 * fx = 0, TF1 * fy = 0, TF2 * fxy = 0);
  ~PDFGen();
};

