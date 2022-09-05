#include "PID_Det.h"

PID_Det::PID_Det(IO_TYPE io_type, TFile *file_input, TFile *file_output,
                 double (*fcn)(double *, double *))
    : PID_Def(io_type, file_input, file_output), mfcn(fcn) {
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
                 double (*fcn)(double *, double *))
    : PID_Def(*pid_def), mfcn(fcn) {
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

  mn_binning_pid = n_bin;
  InitializeHistogram();
}

void PID_Det::InitializeHistogram(int n_bin, double* binning) {
  mbinning_pid = binning;

  mn_binning_pid = n_bin;
  InitializeHistogram();
}

void PID_Det::InitializeHistogram() {
  for (int i = 0; i < 4; i++) {
    mh1_mean[i] = new TH1F(
        Form("mh1_pid_%.2f_%.2f_%d", *(mbinning_pid), *(mbinning_pid + 1), i),
        Form("mh1_pid_%.2f_%.2f_%d", *(mbinning_pid), *(mbinning_pid + 1), i),
        mn_binning_pid, mbinning_pid);
  }

  for (int i = 0; i < mn_binning_pid; i++) {
    // histogram name keep two decimal places
    mv_h1_pid.push_back(mh2_raw->ProjectionY(
        Form("mh1_pid_%.2f_%.2f", *(mbinning_pid), *(mbinning_pid + 1)),
        mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i)),
        mh2_raw->GetXaxis()->FindBin(*(mbinning_pid + i + 1))));
  }

  for (int i = 0; i < mn_binning_pid; i++) {
    mv_f1_pid.push_back(new TF1(
        Form("mh1_pid_%.2f_%.2f", *(mbinning_pid), *(mbinning_pid + 1)), mfcn,
        mh2_raw->GetYaxis()->GetXmin(), mh2_raw->GetYaxis()->GetXmin(), 9));
  }
}

void PID_Det::ShowPIDHistogram() {
  TCanvas *c_pid = new TCanvas("c_Det_PID", "c_Det_PID", 1000, 900);
  int temp = sqrt(mn_binning_pid);
  c_pid->Divide(temp + 1, temp);
  for (int i = 0; i < mn_binning_pid; i++) {
    c_pid->cd(i+1);
    mv_h1_pid[i]->Draw();
  }
}

void PID_Det::HistogramFitting() {}