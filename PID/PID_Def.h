#ifndef PID_Def_h
#define PID_Def_h

#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraphErrors.h"

#include "iostream"
#include "string"

using namespace std;

namespace PID_FCN{
  double func1();
};

class PID_Def {
public:
  enum IO_TYPE { READ = 0, WRITE = 1 };

  PID_Def(IO_TYPE io_type, TFile *file_input, TFile *file_output);
  PID_Def(IO_TYPE io_type, TFile *file_input);
  PID_Def(const PID_Def &pid_def);

  virtual ~PID_Def();

  virtual void GetRawHistogram(TString name) = 0;
  virtual void ShowRawHistogram() = 0;
  virtual void InitializeHistogram(int n_bin, double edge_low,double edge_high) = 0;
  virtual void InitializeHistogram(int n_bin, double* binning) = 0;
  virtual void InitializeTree(TString name) = 0;
  virtual void FillHistogram(double x, double weight = 1) = 0;
  virtual void ShowHistogram() = 0;

  virtual void HistogramFitting() = 0;

protected:
  IO_TYPE mio_type = READ;
  TFile *mfile_input = nullptr;
  TFile *mfile_output = nullptr;

private:
};

#endif