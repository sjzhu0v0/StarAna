#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

const double mean_nsigma_proton_kaon[48] = {
    -0.223709, -0.423411, -0.40938,  -0.641582, -0.893207, -1.01876, -1.36107,
    -1.63507,  -2.0409,   -2.38416,  -2.89433,  -3.43383,  -4.07528, -4.872,
    -5.72069,  -6.49383,  -7.82935,  -9.03003,  -10.5571,  -12.25,   -13.9245,
    -15.5758,  0,         0,         0,         0,         -15.4918, -13.7758,
    -12.1718,  -10.4989,  -8.97858,  -7.68735,  -6.3319,   -5.54946, -4.6893,
    -3.95151,  -3.31492,  -2.78048,  -2.29458,  -1.94676,  -1.59923, -1.27578,
    -1.01032,  -0.718669, -0.605766, -0.409078, -0.255236, -0.308276};

const double width_nsigma_proton_kaon[48] = {
    1.2394,  0.965979, 1.24232, 1.12715, 1.11812, 1.05936, 1.09499, 1.06919,
    1.04013, 1.08014,  1.0756,  1.08618, 1.11172, 1.15261, 1.20319, 1.24701,
    1.38075, 1.43139,  1.53626, 1.59682, 1.60676, 1.49783, 0,       0,
    0,       0,        1.31868, 1.57424, 1.53307, 1.51861, 1.41276, 1.39152,
    1.2544,  1.23351,  1.17132, 1.11741, 1.05211, 1.08225, 1.07472, 1.02282,
    1.05488, 1.09844,  1.12336, 1.10018, 1.03072, 1.06307, 1.1751,  1.02465};

const double mean_nsigma_proton_proton[48] = {
    -0.103328,  -0.111078,  -0.165212,   0.0139826,  -0.0123604, -0.0577394,
    -0.0446237, -0.0490328, -0.0459855,  -0.0395812, -0.0452575, -0.0509856,
    -0.0922701, -0.129425,  -0.117496,   0.0916116,  -0.0481887, 0.0478433,
    0.0387183,  0.0117268,  -0.112149,   -0.301711,  0,          0,
    0,          0,          -0.368912,   -0.0749731, 0.0468328,  0.123085,
    0.194636,   0.105443,   0.245778,    0.0676463,  0.0263151,  0.050435,
    0.0662597,  0.0789152,  0.0731157,   0.0253937,  0.119285,   0.0987428,
    0.0240438,  0.0260706,  -0.00782479, 0.00939157, -0.0132464, -0.0381823};

const double width_nsigma_proton_proton[48] = {
    0.923224, 0.935401, 0.966151, 1.02318,  0.982594, 0.961813, 1.00113,
    0.987146, 0.995684, 0.978816, 1.00105,  0.955686, 0.974778, 0.966967,
    0.960032, 0.95,     0.936611, 0.918305, 0.914343, 0.928203, 0.928219,
    0.995972, 0,        0,        0,        0,        0.96516,  0.95712,
    0.942974, 0.946046, 0.938076, 0.929107, 0.94565,  0.970148, 0.971565,
    0.971069, 0.975036, 0.990846, 0.988921, 0.968986, 0.972263, 0.991144,
    0.977485, 0.955969, 0.956858, 0.992198, 0.979279, 1.0132};

const double mean_m2_kaon[48] = {
    0.1,      0.190059, 0.128647, 0.170057, 0.187141, 0.197349, 0.199445,
    0.210434, 0.213381, 0.21896,  0.220854, 0.224826, 0.227764, 0.231218,
    0.232195, 0.235684, 0.237577, 0.23958,  0.241871, 0.243714, 0,
    0,        0,        0,        0,        0,        0.1,      0.245918,
    0.244243, 0.243117, 0.241175, 0.23927,  0.237085, 0.233922, 0.232309,
    0.228753, 0.227474, 0.222968, 0.22008,  0.21557,  0.209276, 0.204705,
    0.195118, 0.197599, 0.176666, 0.184083, 0.151428, 0.1};

const double width_m2_kaon[48] = {
    0.159535,  0.0631767, 0.266822,  0.0840135,  0.0801955,  0.0697678,
    0.0567242, 0.0554924, 0.0515661, 0.0419031,  0.0371528,  0.0329176,
    0.0281022, 0.0225141, 0.0188038, 0.0157439,  0.0134262,  0.011807,
    0.0118473, 0.0093113, 0,         0,          0,          0,
    0,         0,         0.0125321, 0.00883589, 0.00953622, 0.0114181,
    0.0118463, 0.0144692, 0.0175597, 0.0193517,  0.0243388,  0.0298879,
    0.0296471, 0.0383602, 0.0449873, 0.0532686,  0.0558499,  0.0507254,
    0.0708683, 0.0701024, 0.0889226, 0.0870102,  0.0847109,  0.209163};

const double mean_m2_proton[48] = {
    0.819386, 0.814974, 0.835217, 0.83721,  0.834062, 0.840912, 0.846223,
    0.848219, 0.851643, 0.855742, 0.859028, 0.861079, 0.863881, 0.866409,
    0.868567, 0.871014, 0.873223, 0.875275, 0.877512, 0.879984, 0.88739,
    0,        0,        0,        0,        0,        0,        0.887618,
    0.885346, 0.882992, 0.880269, 0.879198, 0.877531, 0.875475, 0.873368,
    0.86959,  0.867406, 0.864059, 0.862367, 0.859225, 0.855394, 0.850344,
    0.849425, 0.842389, 0.835377, 0.834463, 0.830695, 0.818862};

const double width_m2_proton[48] = {
    0.112479,  0.105804,  0.0992996, 0.0990598, 0.0904534, 0.0804664, 0.0761259,
    0.071586,  0.0638724, 0.0581887, 0.052492,  0.0491002, 0.0451004, 0.040432,
    0.0378404, 0.0346117, 0.0320296, 0.0312151, 0.0315147, 0.0341335, 0.0461715,
    0,         0,         0,         0,         0,         0,         0.0412178,
    0.0346896, 0.0324912, 0.031229,  0.0329267, 0.0343423, 0.0373377, 0.0415004,
    0.0446935, 0.04978,   0.0530329, 0.0582786, 0.0654914, 0.0732514, 0.078164,
    0.0835332, 0.0880609, 0.0905826, 0.0936618, 0.10266,   0.097147};

// const double mean_m2_proton[48] = {};
// const double mean_m2_kaon[48] = {};

bool is_proton_tpc_only(double pq, double nsigma_proton, double cut_n_nsigma) {
  const double edge_low_pq = -2.4;
  const double edge_high_pq = 2.4;
  int i_pq = (pq - edge_low_pq) / (edge_high_pq - edge_low_pq) * 48;
  if (i_pq < 0 || i_pq > 47)
    return false;
  if (nsigma_proton > mean_nsigma_proton_proton[i_pq] +
                          cut_n_nsigma * width_nsigma_proton_proton[i_pq] ||
      nsigma_proton < mean_nsigma_proton_proton[i_pq] -
                          cut_n_nsigma * width_nsigma_proton_proton[i_pq])
    return false;
  return true;
}

TH1F *GetPurity(PID_Det *pid_det, double *pq_binning, int n_binning) {
  TH1F *h1 = new TH1F("h1_purity", "h1_purity", n_binning, pq_binning);
  for (int i = 0; i < n_binning; i++) {
    double pq = h1->GetBinCenter(i + 1);
    double purity = PID_FCN::purity_gaussian(
        -0.2, 2.3, pid_det->GetFittingResult(3)->GetBinContent(i + 1),
        pid_det->GetFittingResult(4)->GetBinContent(i + 1),
        pid_det->GetFittingResult(5)->GetBinContent(i + 1),
        pid_det->GetFittingResult(0)->GetBinContent(i + 1),
        pid_det->GetFittingResult(1)->GetBinContent(i + 1),
        pid_det->GetFittingResult(2)->GetBinContent(i + 1));
    if (purity < 1 && purity > 0)
      h1->SetBinContent(i + 1, purity);
    else
      h1->SetBinContent(i + 1, 0);
  }
  return h1;
}

// #define TPC_PROTON_KAON
#define PROTON_TOF_APPLIED
// #define TOF_KAON

#ifdef PROTON_TOF_APPLIED
void func_fitting(TH1D *h1, TF1 *f1) {
  string name = f1->GetName();
  // get the number behind the last "_" in the name of the function
  int num = atoi(name.substr(name.find_last_of("_") + 1).c_str());
  num--;
  double pars[6] = {100, mean_m2_kaon[num],   width_m2_kaon[num],
                    100, mean_m2_proton[num], width_m2_proton[num]};
  f1->SetParameters(pars);
  f1->SetParLimits(0, 1, 1e8);
  f1->SetParLimits(1, mean_m2_kaon[num] - 5e-3, mean_m2_kaon[num] + 5e-3);
  f1->SetParLimits(2, width_m2_kaon[num] - 5e-3, width_m2_kaon[num] + 5e-3);
  f1->SetParLimits(3, 1, 1e8);
  f1->SetParLimits(4, mean_m2_proton[num] - 5e-3, mean_m2_proton[num] + 5e-3);
  f1->SetParLimits(5, width_m2_proton[num] - 5e-3, width_m2_proton[num] + 5e-3);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma1 = f1->GetParameter(2);
  double mean1 = f1->GetParameter(1);
  double sigma2 = f1->GetParameter(5);
  double mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "", mean1 - sigma1, mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "", mean1 - sigma1, mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "", mean1 - sigma1, mean2 + 2 * sigma2);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (is_proton_tpc_only(TVector3(mom).Mag() * charge, vars_pid[0], 3)) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[3]);
  }
}
#endif

#ifdef TPC_PROTON_KAON
void func_fitting(TH1D *h1, TF1 *f1) {
  string name = f1->GetName();
  // get the number behind the last "_" in the name of the function
  int num = atoi(name.substr(name.find_last_of("_") + 1).c_str());
  num--;
  double pars[6] = {
      100, mean_nsigma_proton_kaon[num],   width_nsigma_proton_kaon[num],
      100, mean_nsigma_proton_proton[num], width_nsigma_proton_proton[num]};
  f1->SetParameters(pars);
  f1->SetParLimits(0, 1, 1e8);
  f1->SetParLimits(1, mean_nsigma_proton_kaon[num] - 5e-2,
                   mean_nsigma_proton_kaon[num] + 5e-2);
  f1->SetParLimits(2, width_nsigma_proton_kaon[num] - 5e-1,
                   width_nsigma_proton_kaon[num] + 5e-1);
  f1->SetParLimits(3, 1, 1e8);
  f1->SetParLimits(4, mean_nsigma_proton_proton[num] - 5e-2,
                   mean_nsigma_proton_proton[num] + 5e-2);
  f1->SetParLimits(5, width_nsigma_proton_proton[num] - 5e-1,
                   width_nsigma_proton_proton[num] + 5e-1);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma1 = f1->GetParameter(2);
  double mean1 = f1->GetParameter(1);
  double sigma2 = f1->GetParameter(5);
  double mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "",
          mean1 + 0.5 * sigma1 < mean2 - 2 * sigma2 ? mean1 + 0.5 * sigma1
                                                    : mean2 - 2 * sigma2,
          mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "",
          mean1 + 0.5 * sigma1 < mean2 - 2 * sigma2 ? mean1 + 0.5 * sigma1
                                                    : mean2 - 2 * sigma2,
          mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "",
          mean1 + 0.5 * sigma1 < mean2 - 2 * sigma2 ? mean1 + 0.5 * sigma1
                                                    : mean2 - 2 * sigma2,
          mean2 + 2 * sigma2);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  double pq = charge * TVector3(mom).Mag();
  h2->Fill(pq, vars_pid[0]);
}

void hist_format(TH1F *h1) {
  h1->GetXaxis()->SetRangeUser(-2.2, 2.2);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.35), 0);
  h1->SetBinContent(h1->GetXaxis()->FindBin(-0.25), 0);
}
#endif

#ifdef TOF_KAON
void func_fitting(TH1D *h1, TF1 *f1) {
  string name = f1->GetName();
  // get the number behind the last "_" in the name of the function
  double pars[6] = {100, 0., 3.73371e-02, 100, 2.20104e-01, 4.83980e-02};
  f1->SetParameters(pars);
  f1->SetParLimits(0, 1, 1e8);
  f1->SetParLimits(1, -5e-2, 5e-2);
  f1->SetParLimits(2, 0., 3.73371e-01);
  f1->SetParLimits(3, 1, 1e8);
  f1->SetParLimits(4, 0.2 - 5e-2, 0.2 + 5e-2);
  f1->SetParLimits(5, 0., 3.73371e-01);
  h1->Fit(f1, "", "", 0 - 0.1, 0 + 0.1);
  double sigma1 = f1->GetParameter(2);
  double mean1 = f1->GetParameter(1);
  double sigma2 = f1->GetParameter(5);
  double mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "",
          mean1 + 0.5 * sigma1 < mean2 - 2 * sigma2 ? mean1 + 0.5 * sigma1
                                                    : mean2 - 2 * sigma2,
          mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "",
          mean1 + 0.5 * sigma1 < mean2 - 2 * sigma2 ? mean1 + 0.5 * sigma1
                                                    : mean2 - 2 * sigma2,
          mean2 + 2 * sigma2);
  sigma1 = f1->GetParameter(2);
  mean1 = f1->GetParameter(1);
  sigma2 = f1->GetParameter(5);
  mean2 = f1->GetParameter(4);
  h1->Fit(f1, "", "", mean1 - 2 * sigma1, mean2 + 2 * sigma2);
}

void filling_h2_pid(TH2F *h2, double *mom, double charge, double *vars_pid) {
  if (vars_pid[1] < 1 && vars_pid[1] > -1) {
    double pq = charge * TVector3(mom).Mag();
    h2->Fill(pq, vars_pid[3]);
  }
}
#endif