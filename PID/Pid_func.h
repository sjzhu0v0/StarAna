#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

// #define TOF_PID_PROTON
// #define TOF_PID_KAON

// #define TPC_PROTON
#define TPC_KAON

#ifdef TOF_PID_PROTON
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
  if (h1_temp->Integral() < 500)
    return;

  double val_center_bin_max = h1_temp->GetBinCenter(h1_temp->GetMaximumBin());
  delete h1_temp;

  double pars[3] = {100, 0.81, 0.12};
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
  // if (vars_pid[0] < 3 && vars_pid[0] > -3) {
  double pq = charge * TVector3(mom).Mag();
  h2->Fill(pq, vars_pid[3]);
  // }
}

void hist_format(TH1F *h1) { h1->GetXaxis()->SetRangeUser(-2.6, 2.6); }
#endif

#ifdef TOF_PID_KAON
void func_fitting(TH1D *h1, TF1 *f1) {
  TH1D *h1_temp = new TH1D();
  h1_temp = (TH1D *)h1->Clone("h1_temp");
  h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(0.1),
                                h1_temp->GetXaxis()->FindBin(0.3));
  double val_center_bin_max = h1_temp->GetBinCenter(h1_temp->GetMaximumBin());
  if (h1_temp->Integral() < 1000) {
    delete h1_temp;
    return;
  }

  double pars[3] = {100, 0.2, 0.12};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 10.);
  f1->SetParLimits(1, 0.1, 0.3);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", val_center_bin_max - 0.1, val_center_bin_max + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (vars_pid[1] < 1 && vars_pid[1] > -1) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[3]);
  }
}

void hist_format(TH1F *h1) {
  h1->GetXaxis()->SetRangeUser(-2.2, 2.2);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.35), 0);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.25), 0);
}
#endif

#ifdef TPC_PROTON
void func_fitting(TH1D *h1, TF1 *f1) {
  // TH1D *h1_temp = new TH1D();
  // h1_temp = (TH1D *)h1->Clone("h1_temp");
  // h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(-1),
  //                               h1_temp->GetXaxis()->FindBin(1));
  // double val_center_bin_max =
  // h1_temp->GetBinCenter(h1_temp->GetMaximumBin()); if (h1_temp->Integral() <
  // 50) {
  //   delete h1_temp;
  //   return;
  // }

  double pars[3] = {100, 0., 1.};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 10.);
  f1->SetParLimits(1, -1, 1);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 2 * sigma, mean + 2 * sigma);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (vars_pid[3] < 0.82 + 0.2 && vars_pid[3] > 0.82 - 0.1) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[0]);
  }
}

void hist_format(TH1F *h1) {
  h1->GetXaxis()->SetRangeUser(-2.2, 2.2);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.35), 0);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.25), 0);
}
#endif

#ifdef TPC_KAON
void func_fitting(TH1D *h1, TF1 *f1) {
  // TH1D *h1_temp = new TH1D();
  // h1_temp = (TH1D *)h1->Clone("h1_temp");
  // h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(-1),
  //                               h1_temp->GetXaxis()->FindBin(1));
  // double val_center_bin_max =
  // h1_temp->GetBinCenter(h1_temp->GetMaximumBin()); if (h1_temp->Integral() <
  // 50) {
  //   delete h1_temp;
  //   return;
  // }
  double val_center_bin_max = h1->GetBinCenter(h1->GetMaximumBin());
  double pars[3] = {100, val_center_bin_max, 1.};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 100.);
  f1->SetParLimits(1, -1 + val_center_bin_max, 1 + val_center_bin_max);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
}

void func_fitting_32(TH1D *h1, TF1 *f1) {
  // TH1D *h1_temp = new TH1D();
  // h1_temp = (TH1D *)h1->Clone("h1_temp");
  // h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(-1),
  //                               h1_temp->GetXaxis()->FindBin(1));
  // double val_center_bin_max =
  // h1_temp->GetBinCenter(h1_temp->GetMaximumBin()); if (h1_temp->Integral() <
  // 50) {
  //   delete h1_temp;
  //   return;
  // }
  double pars[3] = {383., -13.74, 1.6};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 100.);
  f1->SetParLimits(1, -15, -10);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
}

void func_fitting_36(TH1D *h1, TF1 *f1) {
  // TH1D *h1_temp = new TH1D();
  // h1_temp = (TH1D *)h1->Clone("h1_temp");
  // h1_temp->GetXaxis()->SetRange(h1_temp->GetXaxis()->FindBin(-1),
  //                               h1_temp->GetXaxis()->FindBin(1));
  // double val_center_bin_max =
  // h1_temp->GetBinCenter(h1_temp->GetMaximumBin()); if (h1_temp->Integral() <
  // 50) {
  //   delete h1_temp;
  //   return;
  // }
  double pars[3] = {204.7, -8.2, 1.6};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0., 100.);
  f1->SetParLimits(1, -1 - 8.2, 1 - 8.2);
  f1->SetParLimits(0, 1, 1e8);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "", "", mean - 1.5 * sigma, mean + 1.5 * sigma);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (vars_pid[3] > 0.2 && vars_pid[3] < 0.4) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[0]);
  }
}

void hist_format(TH1F *h1) {
  h1->GetXaxis()->SetRangeUser(-2.4, 2.4);
  // h1->SetBinContent(h1->GetXaxis()->FindBin(-0.35), 0);
  // h1->SetBinContent(h1->GetXaxis()->FindBin(-0.25), 0);
}
#endif