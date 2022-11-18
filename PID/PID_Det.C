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

void PID_Det::InitializeTree(TString tree_name, TString p_var_name,
                             TString name_charge, TString names_var[],
                             int n_var) {
  mNvars = n_var;
  mPid_var = new double[mNvars];

  if (mio_type == WRITE) {
    mTree_input = (TTree *)mfile_input->Get(tree_name);
    mTree_input->SetBranchAddress(p_var_name, mP3);
    mTree_input->SetBranchAddress(name_charge, &mCharge);

    for (int i = 0; i < mNvars; i++) {
      mTree_input->SetBranchAddress(names_var[i], (mPid_var + i));
    }
  }
}

void PID_Det::InitializeHistogram(TH1D *hist_raw, int n_var, ...) {
  va_list ap;
  va_start(ap, n_var);
  for (int i = 0; i < n_var / 2; i++) {
    int n_bin = va_arg(ap, int);
    double *binning = va_arg(ap, double *);
    mvec_nbinning.push_back(n_bin);
    mvec_binning.push_back(binning);
  }
  va_end(ap);

  int n_hist = 1;
  for (int i : mvec_nbinning) {
    n_hist *= i;
  }

#ifdef DEBUG
  cout << n_hist << " histograms will be defined. " << endl;
  cout << mvec_nbinning.size() << endl;
#endif

  for (int i = 0; i < n_hist; i++) {
    string name_hist;
    int index_rem = i;
    for (int j = 0; j < mvec_nbinning.size(); j++) {
      name_hist =
          "_" +
          to_string(index_rem % mvec_nbinning[mvec_nbinning.size() - j - 1]) +
          name_hist;
      index_rem = index_rem / mvec_nbinning[mvec_nbinning.size() - j - 1];
    }
    name_hist = "h1_raw" + name_hist;
#ifdef DEBUG
    cout << name_hist << endl;
#endif
    mVh1_pid.push_back((TH1D *)hist_raw->Clone(name_hist.c_str()));
    name_hist = "f1_" + name_hist;
    mVf1_pid.push_back(new TF1(name_hist.c_str(), mfcn,
                               hist_raw->GetXaxis()->GetXmin(),
                               hist_raw->GetXaxis()->GetXmax(), mNpars));
  }
}

void PID_Det::HistogramFilling(
    void HistFilling(TH1D *h1, double *mom, double charge, double *vars_pid),
    void TagCal(double *mom, double charge, double *vars_tag),
    bool CutCal(double *mom, double charge, double *vars_pid)) {
  for (long iEntry = 0; iEntry < mTree_input->GetEntries(); iEntry++) {
    mTree_input->GetEntry(iEntry);
    if (CutCal(mP3, mCharge, mPid_var)) {
      double *vars_tag = new double[mvec_nbinning.size()];
      TagCal(mP3, mCharge, vars_tag);
      int i_hist = 0;
      bool doFill = true;
      for (int i = 0; i < mvec_nbinning[i]; i++) {
        double var_tag = *(vars_tag + i);
        if (var_tag < *(mvec_binning[i]) ||
            var_tag > *(mvec_binning[i] + mvec_nbinning[i])) {
          doFill = false;
          break;
        }
        for (int j = 0; j < mvec_nbinning[i]; j++) {
          if (var_tag > *(mvec_binning[i] + j) &&
              var_tag < *(mvec_binning[i] + j + 1)) {
            i_hist *= mvec_nbinning[i];
            i_hist += j;
          }
        }
      }
      if (doFill) {
        HistFilling(mVh1_pid[i_hist], mP3, mCharge, mPid_var);
      }
    }
  }
}

void PID_Det::ShowHistogram(int i_hist) {
  gStyle->SetOptStat(0);
  if (i_hist == -1) {
    if (mVh1_pid.size() < 100) {
      TCanvas *c1 = new TCanvas("c_RawHist", "c_RawHist", 1000, 1000);
      int n_division = sqrt(mVh1_pid.size()) + 1;
      c1->Divide(n_division, n_division);
      for (int i = 0; i < mVh1_pid.size(); i++) {
        c1->cd(i + 1);
        mVh1_pid[i]->Draw();
        // Draw the index of the histogram at the middle of the pad
        TLatex *latex = new TLatex();
        latex->DrawLatexNDC(0.5, 0.5, to_string(i).c_str());
      }
    } else {
      // Draw all histograms 100 in per canvas
      int n_division = 4;
      int n_canvas = mVh1_pid.size() / pow(n_division, 2);
      if (mVh1_pid.size() % (int)pow(n_division, 2) != 0)
        n_canvas++;

      for (int i = 0; i < n_canvas; i++) {
        TCanvas *c1 = new TCanvas(Form("c_RawHist_%d", i),
                                  Form("c_RawHist_%d", i), 1000, 1000);
        c1->Divide(4, 4);
        for (int j = 0; j < pow(n_division, 2); j++) {
          if (i * pow(n_division, 2) + j < mVh1_pid.size()) {
            c1->cd(j + 1);
            mVh1_pid[i * pow(n_division, 2) + j]->Draw();
            // Draw the index of the histogram at the middle of the pad
            TLatex *latex = new TLatex();
            latex->SetTextSize(0.1);
            latex->DrawLatexNDC(
                0.1, 0.8,
                Form("%.0f:%s", i * pow(n_division, 2) + j,
                     mVh1_pid[i * pow(n_division, 2) + j]->GetName()));
          }
        }
      }
    }
  } else {
    TCanvas *c1 = new TCanvas("c_RawHist", "c_RawHist", 1000, 1000);
    mVh1_pid[i_hist]->Draw();
    c1->Update();
    c1->WaitPrimitive();
  }
}

void PID_Det::HistogramFitting(Option_t *option, Option_t *goption,
                               Double_t xmin, Double_t xmax, int i_which_hist) {
  if (i_which_hist == -1) {
    for (int i = 0; i < mVh1_pid.size(); i++) {
      if (mVh1_pid[i]->Integral() < 50)
        continue;
      cout << "================================================" << endl;
      cout << i << " function fitting starts." << endl;
      cout << "================================================" << endl;
      mVh1_pid[i]->Fit(mVf1_pid[i], option, goption, xmin, xmax);
    }
  } else {
    mVh1_pid[i_which_hist]->Fit(mVf1_pid[i_which_hist], option, goption, xmin,
                                xmax);
  }
}

void PID_Det::HistogramFitting(void FuncFitting(TH1D *, TF1 *),
                               int i_which_hist) {
  if (i_which_hist == -1) {
    for (int i = 0; i < mVh1_pid.size(); i++) {
      if (mVh1_pid[i]->Integral() < 50)
        continue;
      cout << "================================================" << endl;
      cout << i << " function fitting starts." << endl;
      cout << "================================================" << endl;
      FuncFitting(mVh1_pid[i], mVf1_pid[i]);
    }
  } else {
    FuncFitting(mVh1_pid[i_which_hist], mVf1_pid[i_which_hist]);
  }
}

void PID_Det::FittingParInit(int i_which_fcn, double *pars) {
  mVf1_pid[i_which_fcn]->SetParameters(pars);
}

void PID_Det::FittingParInit(double *pars) {
  for (int i = 0; i < mVh1_pid.size(); i++) {
    mVf1_pid[i]->SetParameters(pars);
  }
}

void PID_Det::FittingTuning(void FuncTuning(TF1 *)) {
  for (int i = 0; i < mVh1_pid.size(); i++)
    FittingTuning(i, FuncTuning);
}

void PID_Det::FittingTuning(int i_which_func, void FuncTuning(TF1 *)) {
  FuncTuning(mVf1_pid[i_which_func]);
}

TH1 *PID_Det::FittingResultRequest(
    TH1 *FuncResult(vector<TF1 *> vec_f1, vector<TH1D *> vec_h1,
                    vector<int> mvec_nbinning, vector<double *> mvec_binning)) {
  return FuncResult(mVf1_pid, mVh1_pid, mvec_nbinning, mvec_binning);
}

void PID_Det::SavingResults() {
  mfile_output->cd();
  for (auto ihist : mVh1_pid)
    ihist->Write();
  mfile_output->Close();
}