#include "PID_Det.h"
#include "TSystem.h"

void run(){
  gSystem->Load("build/libPID_Def.a"); 

  TFile* file_input = 
  PID_Det pid_det(PID_Def::WRITE, );
}