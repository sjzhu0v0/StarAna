#ifndef PID_Det_h
#define PID_Det_h

#include "PID_Def.h"

class PID_Det : public PID_Def {
public:
  PID_Det(IO_TYPE io_type, TFile *file_input, TFile *file_output,
          double (*fcn)(const double *, const double *));
  PID_Det(IO_TYPE io_type, TFile *file_input);
  PID_Det(IO_TYPE io_type, const PID_Def *pid_def);
  PID_Det(IO_TYPE io_type, const PID_Def *pid_def,
          double (*fcn)(const double *, const double *));
  ~PID_Det();

  void GetRawHistogram(TString name);
  void ShowRawHistogram();
  void InitializeHistogram(int n_bin, double edge_low, double edge_high);
  void InitializeHistogram(int n_bin, double binning[]);
  void HistogramFit();
  void ShowPIDHistogram();
  void HistogramFitting();

private:
  TString mname = "PID_Det";
  TH2F *mh2_raw = nullptr;
  vector<TH1D*> mv_h1_pid;
  vector<TF1*> mv_f1_pid = {nullptr};
  double* mbinning_pid = nullptr;
  int mn_binning_pid = 0;

  TH1F *mh1_mean[4] = {nullptr}; // e pi k p
  TH1F *mh1_width[4] = {nullptr}; // e pi k p
  
  double (*mfcn)(const double *, const double *) = nullptr;

  void InitializeHistogram();
  void FillHistogram(double x, double weight = 1);
};

#endif