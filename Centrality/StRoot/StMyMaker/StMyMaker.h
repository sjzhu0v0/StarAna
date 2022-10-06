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
  const Double_t getBTofBeta(const UChar_t iCase = 0) const;
  const Bool_t isGoodTrigger() const;
  const Bool_t isGoodEvent() const;
  const Bool_t isGoodTrack() const;
  const Bool_t isElectron() const;
  const Bool_t isPion() const;
  const Bool_t isKaon() const;
  const Bool_t isProton() const;
  const double sDcaxy_cal(TVector3 p, TVector3 dca);

  const TString mInName;
  const TString mOutName;
  TFile *mOutFile;

  StPicoDstMaker *const mPicoDstMaker;
  const StPicoDst *mPicoDst;
  const StPicoEvent *mEvent;
  const StPicoTrack *mTrack;
  const StPicoBTofPidTraits *mBTofPidTraits;
#ifdef MINI_TREE
  TTree *mOutTree;
  Double_t mVz;
  Double_t mVzVpd;
  Int_t mNBtofMatch;
  Int_t mRefMult;
  Int_t mRefMult3;
#endif
  // event
  TH1D *hnevents;

  TH1F *hvz;

  TH2F *hvzvzvpd;
  TH2F *hvzvzvpd_cut;
  TH1F *hbtofmatchmult;

  TH2F *h2NtofRefmult;
  TH2F *h2NtofRefmult_vpdcut;

  TH1F *hrefmult;
  TH1F *hrefmult_vpdcut;
  TH1F *hrefmult_Ntofcut;
  TH1F *hrefmult_Ntofcut_vpdcut;

  TH1F *hrefmult3;
  TH1F *hrefmult3_vpdcut;
  TH1F *hrefmult3_Ntofcut;
  TH1F *hrefmult3_Ntofcut_vpdcut;

  ClassDef(StMyMaker, 1)
};

#endif