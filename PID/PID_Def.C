#include "PID_Def.h"

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
    cerr << "WRITE mode: please use another function to define PID_Def." << endl;
  }
}

//copy PID_Def
PID_Def::PID_Def(const PID_Def &pid_def) {
  mio_type = pid_def.mio_type;
  mfile_input = pid_def.mfile_input;
  mfile_output = pid_def.mfile_output;
}

PID_Def::~PID_Def() {}

