#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

#define TOF_PID

#ifdef TOF_PID
void func_tuning(TF1 *func) {
  double pars[3] = {200, 0.88, 0.14};
  func->SetParameters(pars);
  func->SetParLimits(1, 0.83, 0.93);
  func->SetParLimits(0, 1, 1e7);
}

void func_fitting(TH1D *h1, TF1 *f1) {
  TH1D *h1_temp = new TH1D();
  h1_temp = (TH1D *)h1->Clone("h1_temp");
  h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(0.7),
                                h1_temp->GetXaxis()->GetNbins());
  double val_center_bin_max = h1_temp->GetBinCenter(h1_temp->GetMaximumBin());
  delete h1_temp;

  double pars[3] = {200, 0.8, 0.16};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 10.);
  f1->SetParLimits(1, 0.7, 1.06);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", val_center_bin_max - 0.1, val_center_bin_max + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
}

void func_fitting_56(TH1D *h1, TF1 *f1) {
  TH1D *h1_temp = new TH1D();
  h1_temp = (TH1D *)h1->Clone("h1_temp");
  h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(0.7),
                                h1_temp->GetXaxis()->GetNbins());
  double val_center_bin_max = h1_temp->GetBinCenter(h1_temp->GetMaximumBin());
  delete h1_temp;

  double pars[3] = {7.6, 0.8, 0.14};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 10.);
  f1->SetParLimits(1, 0.7, 1.06);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", 0.48, 1.25);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (vars_pid[0] < 3 && vars_pid[0] > -3) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[3]);
  }
}
#endif

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
#ifdef TOF_PID
  TH2F *h2_template =
      new TH2F("h2_pid_template", "h2_pid_template", n_binning_pq, edge_low_pq,
               edge_high_pq, 250, -0.2, 2.3);
  pid_det->GetRawHistogram(h2_template, "TrackInfo", "Mom", "Charge",
                           names_var_pid, 4, filling_h2_pid);
  pid_det->ShowRawHistogram();
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->HistogramFitting(func_fitting);
  pid_det->HistogramFitting(func_fitting_56,56);
  pid_det->ShowPIDHistogram();
#endif
}

int main() {
  run();
  return 0;
}