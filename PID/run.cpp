#include "Pid_func.h"

void run() {
  R__LOAD_LIBRARY(./build/libPID_Def.so);

  TFile *file_input = new TFile("/data/work/STAR/pAu200_netproton_fluctuation/"
                                "data/result_with_extend_tree.root");
  TFile *file_output = new TFile("./test.root", "recreate");
  PID_Det *pid_det =
      new PID_Det(PID_Def::WRITE, file_input, file_output, PID_FCN::gaus, 3);

  const int n_binning_pq = 60;
  const double edge_low_pq = -3;
  const double edge_high_pq = 3;
  double pq_binning[n_binning_pq + 1] = {0};

  for (int i = 0; i <= n_binning_pq; i++) {
    pq_binning[i] =
        edge_low_pq + (edge_high_pq - edge_low_pq) / (double)n_binning_pq * i;
  }

  // pid_det->GetRawHistogram("hpqnsigmaproton");
  // pid_det->GetRawHistogram("hpqbtof1obetare");
  // pid_det->GetRawHistogram("hpqbtof1obeta");
  TString names_var_pid[4] = {"nSigmaProton", "nSigmaKaon", "nSigmaPion",
                              "BTofM2"};
#ifdef TOF_PID_PROTON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 250, -0.2, 2.3);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->HistogramFitting(func_fitting_56, 56);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 3;
  TString list_names[3] = {"height", "mean", "width"};
  int list_index[3] = {0, 1, 2};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  hist_format(pid_det->GetFittingResult(0));
  hist_format(pid_det->GetFittingResult(1));
  hist_format(pid_det->GetFittingResult(2));
  pid_det->ShowFittingResult();
#endif

#ifdef TOF_PID_KAON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 250, -0.2, 2.3);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 3;
  TString list_names[3] = {"height", "mean", "width"};
  int list_index[3] = {0, 1, 2};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  hist_format(pid_det->GetFittingResult(0));
  hist_format(pid_det->GetFittingResult(1));
  hist_format(pid_det->GetFittingResult(2));
  pid_det->ShowFittingResult();
#endif

#ifdef TPC_PROTON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 200, -20, 20);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 3;
  TString list_names[3] = {"height", "mean", "width"};
  int list_index[3] = {0, 1, 2};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  // hist_format(pid_det->GetFittingResult(0));
  // hist_format(pid_det->GetFittingResult(1));
  // hist_format(pid_det->GetFittingResult(2));
  pid_det->ShowFittingResult();
#endif

#ifdef TPC_KAON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 200, -20, 20);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->HistogramFitting(func_fitting,32);
  pid_det->HistogramFitting(func_fitting,36);
  pid_det->ShowRawHistogram();
  // pid_det->ShowPIDHistogram(33,false);
  pid_det->ShowPIDHistogram();

  const int n_var = 3;
  TString list_names[3] = {"height", "mean", "width"};
  int list_index[3] = {0, 1, 2};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  hist_format(pid_det->GetFittingResult(0));
  hist_format(pid_det->GetFittingResult(1));
  hist_format(pid_det->GetFittingResult(2));
  pid_det->ShowFittingResult();
#endif
}

int main() {
  run();
  return 0;
}