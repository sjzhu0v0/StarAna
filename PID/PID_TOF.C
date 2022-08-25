#include "PID_TOF.h"

PID_TOF::PID_TOF(IO_TYPE io_type, TFile *file_input, TFile *file_output,
                 double (*fcn)(const double *, const double *))
    : PID_Def(io_type, file_input, file_output), mfcn(fcn) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_TOF"
         << endl;
}

PID_TOF::PID_TOF(IO_TYPE io_type, TFile *file_input)
    : PID_Def(io_type, file_input), mfcn(nullptr) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_TOF"
         << endl;
}

PID_TOF::PID_TOF(IO_TYPE io_type, const PID_Def *pid_def)
    : PID_Def(*pid_def), mfcn(nullptr) {
  if ((io_type == READ && mfcn != nullptr) ||
      (io_type == WRITE && mfcn == nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_TOF"
         << endl;
}

PID_TOF::PID_TOF(IO_TYPE io_type, const PID_Def *pid_def,
                 double (*fcn)(const double *, const double *))
    : PID_Def(*pid_def), mfcn(fcn) {
  if ((io_type == READ && mfcn == nullptr) ||
      (io_type == WRITE && mfcn != nullptr))
    cerr << "mfcn is ill-defined. please use another function to construct "
            "PID_TOF"
         << endl;
}

PID_TOF::~PID_TOF() {}

void PID_TOF::GetRawHistogram(TString name) {
  if (mio_type == WRITE) {
    mh2_raw = (TH2F *)mfile_input->Get(name);
    if (mh2_raw == nullptr)
      cerr << "mh2_raw is ill-defined." << endl;
  }
}

void PID_TOF::ShowRawHistogram() {
  if (mh2_raw == nullptr) {
    cerr << "There is no h_raw. Nothing has been drawn." << endl;
  } else {
    TCanvas *c_tof_raw = new TCanvas("c_tof_raw", "c_tof_raw", 900, 900);
    c_tof_raw->cd();
    mh2_raw->Draw();
  }
}

void PID_TOF::InitializeHistogram(int n_bin, double edge_low, double edge_high){
  
}
