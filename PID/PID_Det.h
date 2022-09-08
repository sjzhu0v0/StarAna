#ifndef PID_Det_h
#define PID_Det_h

#include "PID_Def.h"

class PID_Det : public PID_Def {
public:
  PID_Det(IO_TYPE io_type, TFile *file_input, TFile *file_output,
          double (*fcn)(double *, double *), int npars);
  PID_Det(IO_TYPE io_type, TFile *file_input);
  PID_Det(IO_TYPE io_type, const PID_Def *pid_def);
  PID_Det(IO_TYPE io_type, const PID_Def *pid_def,
          double (*fcn)(double *, double *), int npars);
  ~PID_Det();

  void GetRawHistogram(TString name);
  void GetRawHistogram(TH2F *h2_template, TString tree_name, TString p_var_name,
                       TString name_charge, TString names_var[], int n_var,
                       void fill(TH2F *, double *, double, double *));
  void ShowRawHistogram();
  void InitializeHistogram(int n_bin, double edge_low, double edge_high);
  void InitializeHistogram(int n_bin, double *binning);
  void ShowPIDHistogram();
  void ShowPIDHistogram(int i_which_hist,bool is_logy);

  void FittingParInit(double *pars);
  void FittingParInit(int i_which_fcn, double *pars);
  void FittingTuning(void FuncTuning(TF1 *));
  void FittingTuning(int i_which_func, void FuncTuning(TF1 *));

  void HistogramFitting(Option_t *option = "", Option_t *goption = "",
                        Double_t xmin = 0, Double_t xmax = 0, int i_which_hist = -1);
  void HistogramFitting(void FuncFitting(TH1D*,TF1*),int i_which_hist = -1);

private:
  TString mName = "PID_Det";
  TH2F *mh2_raw = nullptr;
  vector<TH1D *> mVh1_pid;
  vector<TF1 *> mVf1_pid;

  double *mbinning_pid = nullptr;
  int mNbinning_pid = 0;
  int mNpars;

  TH1F *mh1_mean[4] = {nullptr};  // e pi k p
  TH1F *mh1_width[4] = {nullptr}; // e pi k p

  double (*mfcn)(double *, double *) = nullptr;

  void InitializeHistogram();
  void FillHistogram(double x, double weight = 1);

  double *mPid_var;
  int mNvars = 0;
};

#endif