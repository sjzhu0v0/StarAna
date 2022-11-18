#ifndef PID_Def_h
#define PID_Def_h

#include "TCanvas.h"
#include "Math/ProbFuncMathCore.h"
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
#include "TH3F.h"
#include "TStyle.h"

#include "iostream"
#include "string"

using namespace std;

namespace PID_FCN {
double gaus(double *x, double *p);
double two_gaus(double *x, double *p);
double four_gaus(double *x, double *p);
double three_gaus(double *x, double *p);
double purity_gaussian(double bound_low = -20, double bound_high = 20,
                     double mean1 = 0, double sigma1 = 1, double height1 = 1,
                     double mean2 = 2, double sigma2 = 1, double height2 = 1);
}; // namespace PID_FCN

class PID_Def {
public:
  enum IO_TYPE { READ = 0, WRITE = 1 };

  PID_Def(IO_TYPE io_type, TFile *file_input, TFile *file_output);
  PID_Def(IO_TYPE io_type, TFile *file_input);
  PID_Def(const PID_Def &pid_def);

  virtual ~PID_Def();

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