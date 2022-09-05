#include "PID_Def.h"
#include "PID_Det.h"
#include "TSystem.h"

void run(){
  gSystem->Load("/home/sjzhu/STAR/Code/StarAna/PID/libPID_def.so"); 

  // double x[1] = {0.};
  // double p[12] = {0.};
  // cout << PID_FCN::four_gaus(x,p) << endl;

  TFile* file_input = new TFile("/data/work/STAR/pAu200_netproton_fluctuation/data/QA.root");
  TFile* file_output = new TFile("./test.root","recreate");
  PID_Det* pid_det = new PID_Det(PID_Def::WRITE,file_input,file_output,PID_FCN::four_gaus);

  const int n_binning_pq = 40;
  const double edge_low_pq = -5.0;
  const double edge_high_pq = 5.0;
  double pq_binning[n_binning_pq+1] = {0};

  for (int i =0; i <= n_binning_pq; i++)
  {
    pq_binning[i] = edge_low_pq + (edge_high_pq - edge_low_pq) / (double)n_binning_pq * i;
  }

  pid_det->GetRawHistogram("hpqnsigmaproton");
  pid_det->ShowRawHistogram();
  pid_det->InitializeHistogram(n_binning_pq,pq_binning);
  pid_det->ShowPIDHistogram();
}