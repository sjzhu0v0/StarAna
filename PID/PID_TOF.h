#ifndef PID_TOF_h
#define PID_TOF_h

#include "PID_Def.h"

class PID_TOF : public PID_Def {
public:
  PID_TOF(IO_TYPE io_type, TFile *file_input, TFile *file_output,
          double (*fcn)(const double *, const double *));
  PID_TOF(IO_TYPE io_type, TFile *file_input);
  PID_TOF(IO_TYPE io_type, const PID_Def *pid_def);
  PID_TOF(IO_TYPE io_type, const PID_Def *pid_def,
          double (*fcn)(const double *, const double *));
  ~PID_TOF();

  void GetRawHistogram(TString name);
  void ShowRawHistogram();
  void InitializeHistogram(int n_bin, double edge_low, double edge_high);
  void InitializeHistogram(int n_bin, double *binning);
  void InitializeTree(TString name);
  void FillHistogram(double x, double weight = 1);
  void ShowHistogram();
  void HistogramFitting();

private:
  TH2F *mh2_raw = nullptr;
  TH1F **mh1_pid = nullptr;
  TF1 **mf1_pid = nullptr;
  TTree *mtree_pid_plus = nullptr;
  TTree *mtree_pid_minus = nullptr;
  double (*mfcn)(const double *, const double *) = nullptr;
};

#endif