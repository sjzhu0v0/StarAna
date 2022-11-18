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
  void InitializeTree(TString tree_name, TString p_var_name,
                      TString name_charge, TString names_var[], int n_var);
  void InitializeHistogram(TH1D *hist_raw, int n_var, ...);
  void HistogramFilling(
      void HistFilling(TH1D *h1, double *mom, double charge, double *vars_pid),
      void TagCal(double *mom, double charge, double *vars_tag),
      bool CutCal(double *mom, double charge, double *vars_pid));
  void ShowHistogram(int i_hist = -1);
  // void ShowPIDHistogram(int i_which_hist, bool is_logy);

  void FittingParInit(double *pars);
  void FittingParInit(int i_which_fcn, double *pars);
  void FittingTuning(void FuncTuning(TF1 *));
  void FittingTuning(int i_which_func, void FuncTuning(TF1 *));

  void HistogramFitting(Option_t *option = "", Option_t *goption = "",
                        Double_t xmin = 0, Double_t xmax = 0,
                        int i_which_hist = -1);
  virtual void HistogramFitting(void FuncFitting(TH1D *, TF1 *),
                                int i_which_hist = -1);

  TH1 *FittingResultRequest(TH1 *FuncResult(vector<TF1 *> vec_f1,
                                            vector<TH1D *> vec_h1,
                                            vector<int> mvec_nbinning,
                                            vector<double *> mvec_binning));
  TH1F *GetFittingResult(int i_which_hist);
  // void ShowFittingResult(int i_which_hist = -1);

  void SavingResults();

private:
  TString mName = "PID_Det";

  vector<int> mvec_nbinning;
  vector<double *> mvec_binning;
  vector<TH1D *> mVh1_pid;
  vector<TF1 *> mVf1_pid;

  double *mbinning_pid = nullptr;
  int mNpars;

  double (*mfcn)(double *, double *) = nullptr;

  void FillHistogram(double x, double weight = 1);

  double *mPid_var;
  int mNvars = 0;
};

#endif
