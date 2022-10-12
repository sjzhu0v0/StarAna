#ifndef StMyMaker_h
#define StMyMaker_h

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"

#include "StChain/StMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#define MINI_TREE

class StMaker;
class StPicoBTofPidTraits;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;

class StMyMaker : public StMaker {
public:
  StMyMaker(const Char_t *const name, const Char_t *const inName,
            const Char_t *const outName, StPicoDstMaker *const picoDstMaker);
  virtual ~StMyMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *option = "");
  virtual Int_t Finish();

private:
  const Int_t MakeEvent();
  const Int_t MakeTrack(const Int_t it);
  const Bool_t isGoodTrigger() const;
  const Bool_t isGoodEvent() const;
  const Bool_t isGoodTrack() const;
  const Bool_t isProton() const;

  const TString mInName;
  const TString mOutName;
  TFile *mOutFile;
  
  StPicoDstMaker *const mPicoDstMaker;
  const StPicoDst *mPicoDst;
  const StPicoEvent *mEvent;
  const StPicoTrack *mTrack;
#ifdef MINI_TREE
  TTree *mOutTree;
  Int_t mRefMult3;
#endif
  Int_t mRunIndex;
  Int_t mNPTracks;
  Int_t mNNTracks;
  Int_t mNProtonP;
  Int_t mNProtonM;

  // event
  TH1D *hnevents;

  ClassDef(StMyMaker, 1)
};

#endif