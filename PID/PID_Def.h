#ifndef PID_Def_h
#define PID_Def_h

#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "iostream"
#include "string"

using namespace std;

namespace PID_FCN {
double gaus(double *x, double *p);
double four_gaus(double *x, double *p);
double three_gaus(double *x, double *p);
}; // namespace PID_FCN

class PID_Def {
public:
  enum IO_TYPE { READ = 0, WRITE = 1 };

  PID_Def(IO_TYPE io_type, TFile *file_input, TFile *file_output);
  PID_Def(IO_TYPE io_type, TFile *file_input);
  PID_Def(const PID_Def &pid_def);

  virtual ~PID_Def();

  virtual void GetRawHistogram(TString name) = 0;
  virtual void ShowRawHistogram() = 0;
  virtual void InitializeHistogram(int n_bin, double edge_low,
                                   double edge_high) = 0;
  virtual void InitializeHistogram(int n_bin, double *binning) = 0;
  virtual void ShowPIDHistogram() = 0;

  virtual void HistogramFitting(Option_t *option, Option_t *goption,
                                Double_t xmin, Double_t xmax,
                                int i_which_hist) = 0;

protected:
  IO_TYPE mio_type = READ;
  TFile *mfile_input = nullptr;
  TFile *mfile_output = nullptr;
  TTree *mTree_input = nullptr;

  double *mP3 = new double[3];
  short mCharge;

private:
};

#endif