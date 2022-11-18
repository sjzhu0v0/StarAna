#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

// #define TPC_PROTON
#define TPC_KAON

#ifdef TPC_PROTON
void HistFilling(TH1D *h1, double *mom, double charge, double *vars_pid) {
  double pq = charge * TVector3(mom).Perp();
  h1->Fill(vars_pid[0]);
}

void TagCal(double *mom, double charge, double *vars_tag) {
  double ptq = TVector3(mom).Perp() * charge;
  double eta = TVector3(mom).PseudoRapidity();
  *vars_tag = eta;
  *(vars_tag + 1) = ptq;
}

bool CutCal(double *mom, double charge, double *vars_pid) {
  return vars_pid[3] < 0.82 + 0.3 && vars_pid[3] > 0.82 - 0.2;
}

void func_fitting(TH1D *h1, TF1 *f1) {
  double pars[3] = {100, 0., 1.};
  f1->SetParameters(pars);
  f1->SetParLimits(2, 0.8, 1.2);
  f1->SetParLimits(1, -0.9, 0.9);
  f1->SetParLimits(0, 1, 1e10);
  h1->Fit(f1, "L", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  h1->GetXaxis()->SetRangeUser(-3, 3);
}

TH1 *MeanRequest(vector<TF1 *> vec_f1, vector<TH1D *> vec_h1,
                 vector<int> mvec_nbinning, vector<double *> mvec_binning) {
  TH2D *h2_mean = new TH2D("h2_mean", "h2_mean", mvec_nbinning[0],
                           mvec_binning[0], mvec_nbinning[1], mvec_binning[1]);
  for (int i = 0; i < vec_f1.size(); i++) {
    int i_x, i_y;
    i_y = i % mvec_nbinning[1];
    i_x = (i / mvec_nbinning[1]) % mvec_nbinning[0];
    h2_mean->SetBinContent(i_x + 1, i_y + 1, vec_f1[i]->GetParameter(1));
  }
  return h2_mean;
}
#endif

#ifdef TPC_KAON
void HistFilling(TH1D *h1, double *mom, double charge, double *vars_pid) {
  double pq = charge * TVector3(mom).Perp();
  h1->Fill(vars_pid[0]);
}

void TagCal(double *mom, double charge, double *vars_tag) {
  double ptq = TVector3(mom).Perp() * charge;
  double eta = TVector3(mom).PseudoRapidity();
  *vars_tag = eta;
  *(vars_tag + 1) = ptq;
}

bool CutCal(double *mom, double charge, double *vars_pid) {
  return vars_pid[3] > 0.15 && vars_pid[3] < 0.3;
}

void func_fitting(TH1D *h1, TF1 *f1) {
  double pars[3] = {100, 0., 1.};

  // pars[1] = h1->GetBinCenter(h1->GetMaximumBin());
  pars[1] = h1->GetMean();

  f1->SetParameters(pars);
  f1->SetParLimits(0, 1, 1e10);
  f1->SetParLimits(1, -0.9 + pars[1], 0.9 + pars[1]);
  f1->SetParLimits(2, 0.8, .2);
  h1->Fit(f1, "L", "", 0 - 0.1, 0 + 0.1);
  double sigma = f1->GetParameter(2);
  double mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  sigma = f1->GetParameter(2);
  mean = f1->GetParameter(1);
  h1->Fit(f1, "L", "", mean - 2 * sigma, mean + 2 * sigma);
  // h1->GetXaxis()->SetRangeUser(-3, 3);
}

TH1 *MeanRequest(vector<TF1 *> vec_f1, vector<TH1D *> vec_h1,
                 vector<int> mvec_nbinning, vector<double *> mvec_binning) {
  TH2D *h2_mean = new TH2D("h2_mean", "h2_mean", mvec_nbinning[0],
                           mvec_binning[0], mvec_nbinning[1], mvec_binning[1]);
  for (int i = 0; i < vec_f1.size(); i++) {
    int i_x, i_y;
    i_y = i % mvec_nbinning[1];
    i_x = (i / mvec_nbinning[1]) % mvec_nbinning[0];
    h2_mean->SetBinContent(i_x + 1, i_y + 1, vec_f1[i]->GetParameter(1));
  }
  return h2_mean;
}
#endif