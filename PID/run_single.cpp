#include "Pid_func_single.h"

void run_single() {
  R__LOAD_LIBRARY(./build/libPID_Def.so);

  TFile *file_input = new TFile("/data/work/STAR/pAu200_netproton_fluctuation/"
                                "data/result_with_extended_tree.root");
  TFile *file_output = new TFile("./test.root", "recreate");
  PID_Det *pid_det =
      new PID_Det(PID_Def::WRITE, file_input, file_output, PID_FCN::gaus, 3);

  const int n_binning_pq = 40;
  const double edge_low_pq = -2.;
  const double edge_high_pq = 2.;

  const int n_binning_y = 10;
  const double edge_low_y = -1.;
  const double edge_high_y = 1.;

  double pq_binning[n_binning_pq + 1] = {0};

  for (int i = 0; i <= n_binning_pq; i++) {
    pq_binning[i] =
        edge_low_pq + (edge_high_pq - edge_low_pq) / (double)n_binning_pq * i;
  }

  double y_binning[n_binning_y + 1] = {0};
  for (int i = 0; i <= n_binning_y; i++) {
    y_binning[i] =
        edge_low_y + (edge_high_y - edge_low_y) / (double)n_binning_y * i;
  }

  TString names_var_pid[4] = {"nSigmaProton", "nSigmaKaon", "nSigmaPion",
                              "BTofM2"};

#ifdef TPC_PROTON
  TH1D *h1_template =
      new TH1D("h1_pid_template", "h1_pid_template", 200, -30, 10);
  pid_det->InitializeTree("TrackInfo", "Mom", "Charge", names_var_pid, 4);
  pid_det->InitializeHistogram(h1_template, 2 * 2, n_binning_y, y_binning,
                               n_binning_pq, pq_binning);
  pid_det->HistogramFilling(HistFilling, TagCal, CutCal);
  pid_det->HistogramFitting(func_fitting);
  TH2D* h2_mean = (TH2D*)pid_det->FittingResultRequest(MeanRequest);
  // pid_det->ShowHistogram();

  TCanvas* c_mean = new TCanvas("c_mean", "c_mean", 800, 600);
  h2_mean->Draw("colz");
#endif

#ifdef TPC_KAON
  TH1D *h1_template =
      new TH1D("h1_pid_template", "h1_pid_template", 200, -30, 10);
  pid_det->InitializeTree("TrackInfo", "Mom", "Charge", names_var_pid, 4);
  pid_det->InitializeHistogram(h1_template, 2 * 2, n_binning_y, y_binning,
                               n_binning_pq, pq_binning);
  pid_det->HistogramFilling(HistFilling, TagCal, CutCal);
  pid_det->HistogramFitting(func_fitting);
  TH2D* h2_mean = (TH2D*)pid_det->FittingResultRequest(MeanRequest);
  pid_det->ShowHistogram();

  TCanvas* c_mean = new TCanvas("c_mean", "c_mean", 800, 600);
  h2_mean->Draw("colz");
#endif

  // pid_det->SavingResults();
}

int main() {
  run_single();
  return 0;
}