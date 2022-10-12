#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <map>
#include <vector>

#include "StBTofUtil/tofPathLength.hh"
#include "StChain/StMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StThreeVector.hh"
#include "phys_constants.h"

#include "StMyCuts.h"
#include "StMyMaker.h"

ClassImp(StMyMaker)
    //-----------------------------------------------------------------------------
    StMyMaker::StMyMaker(const Char_t *const name, const Char_t *const inName,
                         const Char_t *const outName,
                         StPicoDstMaker *const picoDstMaker)
    : StMaker(name), mInName(inName), mOutName(outName),
      mPicoDstMaker(picoDstMaker) {
  mOutFile = nullptr;
  Clear();
}

//-----------------------------------------------------------------------------
StMyMaker::~StMyMaker() {}

//-----------------------------------------------------------------------------
const double StMyMaker::sDcaxy_cal(TVector3 p, TVector3 dca) {
  p = p.Unit();
  float cosl = p.Perp();
  return -p.Y() / cosl * dca.x() + p.x() / cosl * dca.y();
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Init() {
  cout << "Initializing histograms." << endl;

  mOutFile = new TFile(mOutName, "recreate");
  //	mOutFile->cd();

  const Int_t nRunIndices = StMyCuts::RunIdIndex.size();

  hnevents = new TH1D("hnevents", "", 2, -0.5, 1.5);

#ifdef MINI_TREE
  mOutTree = new TTree("TrackInfo", "TrackInfo");
  //branch: RefMult3, I; NProtonP, I; NProtonM, I
  mOutTree->Branch("RefMult3", &mRefMult3, "RefMult3/I");
  mOutTree->Branch("NProtonP", &mNProtonP, "NProtonP/I");
  mOutTree->Branch("NProtonM", &mNProtonM, "NProtonM/I");
#endif
  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Finish() {
  cout << "Writing histograms." << endl;

  mOutFile->cd();

  hnevents->Write();

#ifdef MINI_TREE
  mOutTree->Write();
  cout << mOutTree->GetEntries() << " events have been recorded." << endl;
#endif
  mOutFile->Close();

  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyMaker::Clear(Option_t *option) {
  mPicoDst = nullptr;
  mEvent = nullptr;
  mTrack = nullptr;
  mBTofPidTraits = nullptr;
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Make() {
  if (!mPicoDstMaker) {
    LOG_WARN << "No PicoDstMaker! Skip!" << endm;
    return kStWarn;
  }

  mPicoDst = (const StPicoDst *)mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    LOG_WARN << "No PicoDst! Skip!" << endm;
    return kStWarn;
  }

  MakeEvent();

  return kStOK;
}

//-----------------------------------------------------------------------------
const Int_t StMyMaker::MakeEvent() {
  mEvent = (const StPicoEvent *)mPicoDst->event();
  if (!mEvent) {
    LOG_WARN << "Error opening picoDst Event! Skip!" << endm;
    return kStWarn;
  }

  hnevents->Fill(0.);

  if (!isGoodTrigger() || !isGoodEvent())
    return kStOK;

  mRunIndex = StMyCuts::RunIdIndex.at(mEvent->runId());
  const TVector3 PrimaryVertex(mEvent->primaryVertex().x(),
                               mEvent->primaryVertex().y(),
                               mEvent->primaryVertex().z());

  hnevents->Fill(1.);

  mRefMult3 = mEvent->refMult3();

  mNPTracks = mNNTracks = mNProtonP = mNProtonM = 0;
  const Int_t nTracks = mPicoDst->numberOfTracks();
  for (Int_t it = 0; it < nTracks; it++)
    MakeTrack(it);

#ifdef MINI_TREE
  mOutTree->Fill();
#endif

  return kStOK;
}

//-----------------------------------------------------------------------------
const Int_t StMyMaker::MakeTrack(const Int_t it) {
  mTrack = (const StPicoTrack *)mPicoDst->track(it);
  if (!mTrack) {
    LOG_WARN << "Error opening picoDst Track! Skip!" << endm;
    return kStWarn;
  }
  mBTofPidTraits = mTrack->bTofPidTraitsIndex() != -1
                       ? (const StPicoBTofPidTraits *)mPicoDst->btofPidTraits(
                             mTrack->bTofPidTraitsIndex())
                       : nullptr;

  if (!isGoodTrack())
    return kStOK;

  const TVector3 Dca = mTrack->gDCA(mEvent->primaryVertex());
  const TVector3 PMom = mTrack->pMom();
  const Double_t SDcaXY = sDcaxy_cal(PMom, Dca);
  const Double_t PPt = PMom.Pt();

  if (mTrack->charge() > 0) {
    mNPTracks++;
    if (isProton())
      mNProtonP++;
  }
  if (mTrack->charge() < 0) {
    mNNTracks++;
    if (isProton())
      mNProtonM++;
  }

  return kStOK;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodTrigger() const {
  for (auto iTrigger : StMyCuts::Trigger)
    if (mEvent->isTrigger(iTrigger))
      return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodEvent() const {
  const Double_t Vx = mEvent->primaryVertex().x();
  const Double_t Vy = mEvent->primaryVertex().y();
  const Double_t Vz = mEvent->primaryVertex().z();
  if (TMath::Abs(Vx) < 1.e-5 && TMath::Abs(Vy) < 1.e-5 &&
      TMath::Abs(Vz) < 1.e-5)
    return kFALSE;
  if (Vz <= StMyCuts::VzCut[0] || Vz >= StMyCuts::VzCut[1])
    return kFALSE;
  // cut:vzvpd - vz
  if (TMath::Abs(mEvent->vzVpd() - Vz) > StMyCuts::cut_VzMinusVzvpd)
    return kFALSE;
  // rejecting bad runs
  for (auto iRun : StMyCuts::list_badrun)
    if (mRunIndex == iRun)
      return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodTrack() const {
  if (!mTrack->isPrimary())
    return kFALSE;
  if (mTrack->gDCA(mEvent->primaryVertex()).Mag() >= StMyCuts::DcaCut)
    return kFALSE;
  if (mTrack->nHitsFit() <= StMyCuts::NHitsFitCut)
    return kFALSE;
  if ((const Double_t)mTrack->nHitsFit() / mTrack->nHitsMax() <=
      StMyCuts::NHitsRatioCut)
    return kFALSE;
  if (mTrack->pMom().Perp() < StMyCuts::PtCutQ[0] ||
      mTrack->pMom().Perp() > StMyCuts::PtCutQ[1])
    return kFALSE;
  if (mTrack->pMom().PseudoRapidity() < StMyCuts::EtaCutQ[0] ||
      mTrack->pMom().PseudoRapidity() > StMyCuts::EtaCutQ[1])
    return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isProton() const {
  if (TMath::Abs(mTrack->nSigmaProton()) >= StMyCuts::NSigmaCut)
    return kFALSE;
  return kTRUE;
}
