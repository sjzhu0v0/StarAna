#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

#define TOF_PID

#ifdef TOF_PID
void func_tuning(TF1* func){
  func->SetParLimits(2,-0.1,0.1);
  func->SetParLimits(5,0.18,0.22);
  func->SetParLimits(7,0.7,0.94);
  func->SetParLimits(0,10,1e8);
  func->SetParLimits(3,10,1e8);
  func->SetParLimits(6,10,1e8);
}
#endif

void run() {
  R__LOAD_LIBRARY(./build/libPID_Def.so);

  TFile *file_input =
      new TFile("/data/work/STAR/pAu200_netproton_fluctuation/data/result_with_minitree.root");
  TFile *file_output = new TFile("./test.root", "recreate");
  PID_Det *pid_det =
      new PID_Det(PID_Def::WRITE, file_input, file_output, PID_FCN::three_gaus,9);

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
#ifdef TOF_PID
  pid_det->GetRawHistogram("hpqbtofm2");
  pid_det->ShowRawHistogram();
  pid_det->InitializeHistogram(n_binning_pq, pq_binning);
  double pars[9] = {1000,0.01,0.01,500,0.25,0.01,200,0.85,0.08};
  pid_det->FittingParInit(pars);
  pid_det->FittingTuning(func_tuning);
  pid_det->HistogramFitting();
  pid_det->ShowPIDHistogram();
#endif

}

int main() {
  run();
  return 0;
}