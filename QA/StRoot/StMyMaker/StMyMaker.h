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

  TTree *mOutTree;
  Double_t mMom_Minitree[3];
  Short_t mCharge_Minitree;
  Double_t mNSigmaProton_Minitree;
  Double_t mBTofM2_Minitree_Minitree;

  Int_t mRunIndex;
  Int_t mNPTracks;
  Int_t mNNTracks;
  Int_t mNElectronP;
  Int_t mNElectronM;
  Int_t mNPionP;
  Int_t mNPionM;
  Int_t mNKaonP;
  Int_t mNKaonM;
  Int_t mNProtonP;
  Int_t mNProtonM;
  Double_t mQ1xTpc;
  Double_t mQ1yTpc;
  Double_t mQ2xTpc;
  Double_t mQ2yTpc;

  // event
  TH1D *hnevents;
  TH2F *hvxvy;
  TH1F *hvz;
  TH2F *hvzvzvpd;
  TH1F *hnvpdhitseast;
  TH1F *hnvpdhitswest;
  TH1F *hzdcx;
  TH1F *hrefmult;
  TH1F *hrefmult2;
  TH1F *hrefmult3;
  TH1F *hrefmult4;
  TH1F *hgrefmult;
  TH1F *hbtofmatchmult;
  TH1F *hbtoftraymult;
  TH1F *hnptracks;
  TH1F *hnntracks;
  TH1F *hnelectronp;
  TH1F *hnelectronm;
  TH1F *hnpionp;
  TH1F *hnpionm;
  TH1F *hnkaonp;
  TH1F *hnkaonm;
  TH1F *hnprotonp;
  TH1F *hnprotonm;

  TProfile *pvx;
  TProfile *pvy;
  TProfile *pvz;
  TProfile *pvr;
  TProfile *pvzmvzvpd;
  TProfile *pnvpdhitseast;
  TProfile *pnvpdhitswest;
  TProfile *pzdcx;
  TProfile *prefmult;
  TProfile *prefmult2;
  TProfile *prefmult3;
  TProfile *prefmult4;
  TProfile *pgrefmult;
  TProfile *pbtofmatchmult;
  TProfile *pbtoftraymult;
  TProfile *pnptracks;
  TProfile *pnntracks;
  TProfile *pnelectronp;
  TProfile *pnelectronm;
  TProfile *pnpionp;
  TProfile *pnpionm;
  TProfile *pnkaonp;
  TProfile *pnkaonm;
  TProfile *pnprotonp;
  TProfile *pnprotonm;
  TProfile *pq1xtpc;
  TProfile *pq1ytpc;
  TProfile *pq2xtpc;
  TProfile *pq2ytpc;

  // track
  TH1F *hdca;
  TH1F *hdcaxy;
  TH1F *hdcaphi;
  TH1F *hdcaz;
  TH1F *hsdcaxy;
  TH1F *hnhitsfit;
  TH1F *hnhitsratio;
  TH1F *hnhitsdedx;
  TH1F *hpt;
  TH1F *hphi;
  TH1F *heta;
  TH2F *hpqdedx;
  TH2F *hpqnsigmaelectron;
  TH2F *hpqnsigmapion;
  TH2F *hpqnsigmakaon;
  TH2F *hpqnsigmaproton;
  TH2F *hbtofyzlocal;
  TH2F *hbtof1obetacheck;
  TH2F *hpqbtof1obeta;
  TH2F *hpqbtof1obetare;
  TH2F *hpqbtofm2;
  TH2F *hpqbtofm2re;

  TProfile *pdca;
  TProfile *pdcaxy;
  TProfile *pdcaphi;
  TProfile *pdcaz;
  TProfile *psdcaxy;
  TProfile *psdcaxy2;
  TProfile *pnhitsfit;
  TProfile *pnhitsratio;
  TProfile *pnhitsdedx;
  TProfile *ppt;
  TProfile *pphi;
  TProfile *peta;
  TProfile *pdedx;
  TProfile *pbtof1obeta;
  TProfile *pbtofylocal;
  TProfile *pbtofzlocal;

  ClassDef(StMyMaker, 1)
};

#endif