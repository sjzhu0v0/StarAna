#include "PID_Def.h"

double PID_FCN::gaus(double *x, double *p) {
  return p[0] * TMath::Gaus(x[0], p[1], p[2]);
}

double PID_FCN::two_gaus(double *x, double *p) {
  return p[0] * TMath::Gaus(x[0], p[1], p[2]) +
         p[3] * TMath::Gaus(x[0], p[4], p[5]);
}

double PID_FCN::three_gaus(double *x, double *p) {
  return p[0] * TMath::Gaus(x[0], p[1], p[2]) +
         p[3] * TMath::Gaus(x[0], p[4], p[5]) +
         p[6] * TMath::Gaus(x[0], p[7], p[8]);
}

double PID_FCN::four_gaus(double *x, double *p) {
  return p[0] * TMath::Gaus(x[0], p[1], p[2]) +
         p[3] * TMath::Gaus(x[0], p[4], p[5]) +
         p[6] * TMath::Gaus(x[0], p[7], p[8]) +
         p[9] * TMath::Gaus(x[0], p[10], p[11]);
}

double PID_FCN::purity_gaussian(double bound_low, double bound_high,
                                double height1, double mean1, double sigma1,
                                double height2, double mean2, double sigma2) {
  // cout input
  cout << "bound_low: " << bound_low << " bound_high: " << bound_high
       << " mean1: " << mean1 << " sigma1: " << sigma1
       << " height1: " << height1 << " mean2: " << mean2
       << " sigma2: " << sigma2 << " height2: " << height2 << endl;

  // two gaussian func with different mean, sigma and height
  TF1 *f1 = new TF1("f1_purity_gaussian", gaus, bound_low, bound_high, 3);
  f1->SetParameters(height1, mean1, sigma1);
  TF1 *f2 = new TF1("f2_purity_gaussian", gaus, bound_low, bound_high, 3);
  f2->SetParameters(height2, mean2, sigma2);

  // calculate the purity of f1 in the 3-sigma region
  double n_imperfection =
      ROOT::Math::gaussian_cdf(mean1 - 3 * sigma1, sigma2, mean2) +
      ROOT::Math::gaussian_cdf(mean1 + 3 * sigma1, sigma2, mean2);
  n_imperfection *= f2->Integral(mean1 - 3 * sigma1, mean1 + 3 * sigma1);
  double purity =
      f1->Integral(mean1 - 3 * sigma1, mean1 + 3 * sigma1) /
      (f1->Integral(mean1 - 3 * sigma1, mean1 + 3 * sigma1) + n_imperfection);
  cout << "purity: " << purity
       << " f1 integral: " << f1->Integral(bound_low, bound_high)
       << " f2 integral: " << f2->Integral(bound_low, bound_high)
       << " n_imperfection: " << n_imperfection << endl;
  delete f1;
  delete f2;
  return purity;
}

PID_Def::PID_Def(IO_TYPE io_type, TFile *file_input, TFile *file_output)
    : mfile_input(file_input), mfile_output(file_output), mio_type(io_type) {
  if (mio_type == WRITE)
    if (mfile_input == nullptr || mfile_output == nullptr)
      cerr << "mfile_input or mfile_output is ill-defined." << endl;

  if (mio_type == READ) {
    cerr << "READ mode: please use another function to define PID_Def." << endl;
  }
}

PID_Def::PID_Def(IO_TYPE io_type, TFile *file_input)
    : mfile_input(file_input), mio_type(io_type) {
  if (mio_type == READ)
    if (mfile_input == nullptr)
      cerr << "mfile_input is ill-defined." << endl;

  if (mio_type == WRITE) {
    cerr << "WRITE mode: please use another function to define PID_Def."
         << endl;
  }
}

// copy PID_Def
PID_Def::PID_Def(const PID_Def &pid_def) {
  mio_type = pid_def.mio_type;
  mfile_input = pid_def.mfile_input;
  mfile_output = pid_def.mfile_output;
}

PID_Def::~PID_Def() {}
