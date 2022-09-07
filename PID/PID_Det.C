#include "PID_Det.h"

// #define DEBUG

PID_Det::PID_Det(IO_TYPE io_type, TFile *file_input, TFile *file_output,
                 double (*fcn)(double *, double *), int npars)
    : PID_Def(io_type, file_input, file_output), mfcn(fcn), mNpars(npars) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_Det"
         << endl;
}

PID_Det::PID_Det(IO_TYPE io_type, TFile *file_input)
    : PID_Def(io_type, file_input), mfcn(nullptr) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_Det"
         << endl;
}

PID_Det::PID_Det(IO_TYPE io_type, const PID_Def *pid_def)
    : PID_Def(*pid_def), mfcn(nullptr) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_Det"
         << endl;
}

PID_Det::PID_Det(IO_TYPE io_type, const PID_Def *pid_def,
                 double (*fcn)(double *, double *), int npars)
    : PID_Def(*pid_def), mfcn(fcn), mNpars(npars) {
  if ((io_type == READ && mfcn == nullptr) ||
      (io_type == WRITE && mfcn != nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_Det"
         << endl;
}

PID_Det::~PID_Det() {}

void PID_Det::GetRawHistogram(TString name) {
  if (mio_type == WRITE) {
    mh2_raw = (TH2F *)mfile_input->Get(name);
    mh2_raw->SetName("mh2_raw");
    if (mh2_raw == nullptr)
      cerr << "mh2_raw is ill-defined." << endl;
  }
}

void PID_Det::GetRawHistogram(const TH2F &h2_template, TString tree_name,
                              TString p_var_name, TString name_charge,
                              TString tpc_var_name, TString tof_var_name,
                              TString cut) {
  if (mio_type == WRITE) {
    mh2_raw = new TH2F(h2_template);
    mh2_raw->SetName("mh2_raw");

    mTree_input = (TTree *)mfile_input->Get(tree_name);
    mTree_input->SetBranchAddress(p_var_name, mP3);
    mTree_input->SetBranchAddress(name_charge, &mCharge);

    if (mh2_raw == nullptr)
      cerr << "mh2_raw is ill-defined." << endl;
  }
}

void PID_Det::ShowRawHistogram() {
  if (mh2_raw == nullptr) {
    cerr << "There is no h_raw. Nothing has been drawn." << endl;
  } else {
    TCanvas *c_Det_raw = new TCanvas("c_Det_raw", "c_Det_raw", 900, 900);
    c_Det_raw->cd();
    gPad->SetLogz();
    mh2_raw->Draw("colz");
  }
}

void PID_Det::InitializeHistogram(int n_bin, double edge_low,
                                  double edge_high) {
  double array_temp[(const int)n_bin + 1];
  for (int i = 0; i < n_bin; i++) {
    array_temp[i] = edge_low + (edge_high - edge_low) / (double)n_bin * i;
  }
  mbinning_pid = array_temp;

  mNbinning_pid = n_bin;
  InitializeHistogram();
}

void PID_Det::InitializeHistogram(int n_bin, double *binning) {
  mbinning_pid = binning;

  mNbinning_pid = n_bin;
  InitializeHistogram();
}

void PID_Det::InitializeHistogram() {
  for (int i = 0; i < 4; i++) {
    mh1_mean[i] = new TH1F(Form("mh1_mean_%d", i), Form("mh1_mean_%d", i),
                           mNbinning_pid, mbinning_pid);
    mh1_width[i] = new TH1F(Form("mh1_width_%d", i), Form("mh1_width_%d", i),
                            mNbinning_pid, mbinning_pid);
  }

  // #ifdef DEBUG
  cout << "mv_binning_pid: " << mNbinning_pid << endl;
  // #endif

  for (int i = 0; i < mNbinning_pid; i++) {
    // histogram name keep two decimal places
    mVh1_pid.push_back(mh2_raw->ProjectionY(
        Form("mh1_pid_%d_%d", i, i + 1),
        mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i)),
        mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i + 1))));
#ifdef DEBUG
    cout << "N bins: " << mh2_raw->GetXaxis()->GetNbins() << endl;
    cout << Form("mh1_pid_%d_%d", i, i + 1) << ":"
         << mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i)) << " to "
         << mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i + 1)) << endl;
#endif
  }

  for (int i = 0; i < mNbinning_pid; i++) {
    mVf1_pid.push_back(new TF1(Form("mf1_pid_%d_%d", i, i + 1), mfcn,
                               mVh1_pid[i]->GetXaxis()->GetXmin(),
                               mVh1_pid[i]->GetXaxis()->GetXmax(), mNpars));
  }
}

void PID_Det::ShowPIDHistogram(int i_which_hist) {
  TCanvas *c_pid = new TCanvas("c_Det_PID", "c_Det_PID", 1000, 900);
  c_pid->cd();
  gPad->SetLogy();
  mVh1_pid[i_which_hist]->Draw();
}

void PID_Det::ShowPIDHistogram() {
  TCanvas *c_pid = new TCanvas("c_Det_PID", "c_Det_PID", 1000, 900);
  int temp = sqrt(mNbinning_pid);
  int n_high, n_width;
  if (temp % 2 == 0) {
    n_width = temp;
    n_high = temp + 1;
  } else {
    n_width = temp + 1;
    n_high = temp;
  }

  c_pid->Divide(n_width, n_high);
  for (int i = 0; i < mNbinning_pid / 2; i++) {
    int i2plotted = i + 1 + (n_width / 2) * (i / (n_width / 2));
#ifdef DEBUG
    cout << i2plotted << " " << mNbinning_pid / 2 - i;
#endif
    c_pid->cd(i2plotted);
    gPad->SetLogy();
    mVh1_pid[mNbinning_pid / 2 - i - 1]->Draw();

    TLatex *latex1 =
        new TLatex(0.15, 0.8,
                   Form("%d:%.2f < pq < %.2f", mNbinning_pid / 2 - i - 1,
                        *(mbinning_pid + mNbinning_pid / 2 - i - 1),
                        *(mbinning_pid + mNbinning_pid / 2 - i)));
    latex1->SetNDC();
    latex1->SetTextSize(0.1);
    latex1->SetTextColor(kRed);
    latex1->Draw();

    i2plotted = i + 1 + (n_width / 2) * (i / (n_width / 2) + 1);
    c_pid->cd(i2plotted);
    gPad->SetLogy();
    mVh1_pid[mNbinning_pid / 2 + i]->Draw();
#ifdef DEBUG
    cout << "  " << i2plotted << "  " << mn_binning_pid / 2 + 1 + i << endl;
#endif
    TLatex *latex2 =
        new TLatex(0.15, 0.8,
                   Form("%d:%.2f < pq < %.2f", mNbinning_pid / 2 + i,
                        *(mbinning_pid + mNbinning_pid / 2 + i),
                        *(mbinning_pid + mNbinning_pid / 2 + i + 1)));
    latex2->SetNDC();
    latex2->SetTextSize(0.1);
    latex2->SetTextColor(kRed);
    latex2->Draw();
  }
}

void PID_Det::HistogramFitting() {
  for (int i = 0; i < mNbinning_pid; i++) {
    if (mVh1_pid[i]->Integral() < 1000)
      continue;
    cout << "================================================" << endl;
    cout << i << " function fitting starts." << endl;
    cout << "================================================" << endl;
    mVh1_pid[i]->Fit(mVf1_pid[i], "L");
  }
}

void PID_Det::FittingParInit(int i_which_fcn, double *pars) {
  mVf1_pid[i_which_fcn]->SetParameters(pars);
}

void PID_Det::FittingParInit(double *pars) {
  for (int i = 0; i < mNbinning_pid; i++) {
    mVf1_pid[i]->SetParameters(pars);
  }
}

void PID_Det::FittingTuning(void FuncTuning(TF1 *)) {
  for (int i = 0; i < mNbinning_pid; i++)
    FittingTuning(i, FuncTuning);
}

void PID_Det::FittingTuning(int i_which_func, void FuncTuning(TF1 *)) {
  FuncTuning(mVf1_pid[i_which_func]);
}