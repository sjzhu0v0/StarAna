#include "Pid_func_multi.h"

void run_multi() {
  R__LOAD_LIBRARY(./build/libPID_Def.so);

  TFile *file_input = new TFile("/data/work/STAR/pAu200_netproton_fluctuation/"
                                "data/result_with_extended_tree.root");
  TFile *file_output = new TFile("./test.root", "recreate");
// #if (defined(TPC_PROTON_KAON) || defined(TOF_KAON))
  PID_Det *pid_det = new PID_Det(PID_Def::WRITE, file_input, file_output,
                                 PID_FCN::two_gaus, 6);
// #else
//   PID_Det *pid_det = new PID_Det(PID_Def::WRITE, file_input, file_output,
//                                  PID_FCN::three_gaus, 6);
// #endif

  const int n_binning_pq = 48;
  const double edge_low_pq = -2.4;
  const double edge_high_pq = 2.4;
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

#ifdef PROTON_TOF_APPLIED
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 250, -0.2, 2.3);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 6;
  TString list_names[n_var] = {"height_kaon",   "mean_kaon",   "width_kaon",
                               "height_proton", "mean_proton", "width_proton"};
  int list_index[6] = {0, 1, 2, 3, 4, 5};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  pid_det->ShowFittingResult();

  TH1F *h1_purity = GetPurity(pid_det, pq_binning, n_binning_pq);

  // draw purity c1_purity
  TCanvas *c1_purity = new TCanvas("c1_purity", "c1_purity", 800, 600);
  c1_purity->cd();
  h1_purity->Draw();
#endif

#ifdef TPC_PROTON_KAON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 200, -30, 10);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 6;
  TString list_names[n_var] = {"height_kaon",   "mean_kaon",   "width_kaon",
                               "height_proton", "mean_proton", "width_proton"};
  int list_index[6] = {0, 1, 2, 3, 4, 5};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  // hist_format(pid_det->GetFittingResult(0));
  // hist_format(pid_det->GetFittingResult(1));
  // hist_format(pid_det->GetFittingResult(2));
  pid_det->ShowFittingResult();

  TH1F *h1_purity = GetPurity(pid_det, pq_binning, n_binning_pq);

  // draw purity c1_purity
  TCanvas *c1_purity = new TCanvas("c1_purity", "c1_purity", 800, 600);
  c1_purity->cd();
  h1_purity->GetXaxis()->SetRangeUser(-1.2, 1.2);
  h1_purity->Draw();
#endif

#ifdef TOF_KAON
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 250, -0.2, 2.3);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->ShowRawHistogram();
  pid_det->ShowPIDHistogram();

  const int n_var = 6;
  TString list_names[n_var] = {"height_kaon",   "mean_kaon",   "width_kaon",
                               "height_proton", "mean_proton", "width_proton"};
  int list_index[6] = {0, 1, 2, 3, 4, 5};
  pid_det->FittingResultRequest(n_var, list_index, list_names);
  pid_det->ShowFittingResult();

  TH1F *h1_purity = GetPurity(pid_det, pq_binning, n_binning_pq);

  // draw purity c1_purity
  TCanvas *c1_purity = new TCanvas("c1_purity", "c1_purity", 800, 600);
  c1_purity->cd();
  h1_purity->Draw();
#endif
}

int main() {
  run_multi();
  return 0;
}
