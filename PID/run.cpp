#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

void run() {
  R__LOAD_LIBRARY(./build/libPID_Def.so);

  TFile *file_input =
      new TFile("/data/work/STAR/pAu200_netproton_fluctuation/data/QA.root");
  TFile *file_output = new TFile("./test.root", "recreate");
  PID_Det *pid_det =
      new PID_Det(PID_Def::WRITE, file_input, file_output, PID_FCN::four_gaus);

  const int n_binning_pq = 30;
  const double edge_low_pq = -1.5;
  const double edge_high_pq = 1.5;
  double pq_binning[n_binning_pq + 1] = {0};

  for (int i = 0; i <= n_binning_pq; i++) {
    pq_binning[i] =
        edge_low_pq + (edge_high_pq - edge_low_pq) / (double)n_binning_pq * i;
  }

  // pid_det->GetRawHistogram("hpqnsigmaproton");
  pid_det->GetRawHistogram("hpqbtofm2");
  // pid_det->GetRawHistogram("hpqbtof1obetare");
  // pid_det->GetRawHistogram("hpqbtof1obeta");
  pid_det->ShowRawHistogram();
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  pid_det->ShowPIDHistogram();
}

int main() {
  run();
  return 0;
}