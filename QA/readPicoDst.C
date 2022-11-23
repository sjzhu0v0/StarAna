#include "TChain.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

#include "StChain/StChain.h"

#include "StMyMaker.cxx"

#include "iostream"
#include "time.h"

using namespace std;

class StChain;
class StPicoDstMaker;
class StMyMaker;

void readPicoDst(const Char_t *const inputFile = "test.list",
                 const Char_t *const outputFile = "test.root",
                 Long64_t nEvents = -1) {
  TStopwatch *const stopWatch = new TStopwatch();
  stopWatch->Start();

  const TString SL_version = "SL18h";
  const TString env_SL = gSystem->Getenv("STAR");
  if (!env_SL.Contains(SL_version)) {
    cout << "Environment Star Library does not match the requested library in "
            "readPicoDst.C. Exiting..."
         << endl;
    gSystem->Exit(1);
  }

  gSystem->Load("StarClassLibrary");
  gSystem->Load("StarRoot");
  gSystem->Load("libStDb_Tables");
  gSystem->Load("libgeometry_Tables");
  gSystem->Load("RTS");
  gSystem->Load("StTableUtilities");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");
  gSystem->Load("StTreeMaker");
  gSystem->Load("StIOMaker");
  gSystem->Load("StTriggerDataMaker");
  gSystem->Load("StBichsel");
  gSystem->Load("StEvent");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StDbLib");
  gSystem->Load("StEmcUtil");
  gSystem->Load("StTofUtil");
  gSystem->Load("StPmdUtil");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StMagF");
  gSystem->Load("StMtdHitMaker");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StBTofUtil");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  gSystem->Load("StMyMaker");

  StChain *const chain = new StChain();
  StPicoDstMaker *const picoDstMaker =
      new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDstMaker");
  StMyMaker *const myMaker =
      new StMyMaker("myMaker", inputFile, outputFile, picoDstMaker);

  chain->Init();
  const Long64_t nEntries = picoDstMaker->chain()->GetEntries();
  if (nEvents < 0 || nEvents > nEntries)
    nEvents = nEntries;
  cout << nEntries << " events in chain, " << nEvents << " will be read."
       << endl;
  for (Long64_t ie = 0; ie < nEvents; ie++) {
    if (ie % 100000 == 0) {
      cout << "Working on eventNumber " << ie << ". Time:";
      time_t t = time(NULL);
      char ch[64] = {0};
      char time_str[100] = {0};
      strftime(ch, sizeof(ch) - 1, "%Y%m%d--%H%M%S", localtime(&t));
      sprintf(time_str, "%s", ch);
      cout << time_str << endl;
    }
    chain->Clear();
    chain->Make();
  }
  // chain->EventLoop(nEvents);
  chain->Finish();
  delete chain;

  stopWatch->Stop();
  stopWatch->Print();
}
