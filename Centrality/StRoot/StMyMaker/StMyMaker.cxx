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

  hnevents = new TH1D("hnevents", "", 2, -0.5, 1.5);

  hvz = new TH1F("hvz", ";#it{V}_{#it{z}} (cm);#it{N}_{events}", 80, -100, 100);

  hvzvzvpd =
      new TH2F("hvzvzvpd", ";#it{V}_{#it{z}} (cm);#it{V}_{#it{z}}^{VPD} (cm)",
               80, -100., 100., 80, -100., 100.);
  hvzvzvpd_cut = new TH2F(
      "hvzvzvpd_cut",
      "|V_{z}-V_{zvpd}|<6 cm;#it{V}_{#it{z}} (cm);#it{V}_{#it{z}}^{VPD} (cm)",
      80, -100., 100., 80, -100., 100.);
  hbtofmatchmult = new TH1F("hbtofmatchmult", ";bTofMatchMult;#it{N}_{events}",
                            80, -0.5, 79.5);
  h2NtofRefmult = new TH2F("h2NtofRefmult", ";bTofMatchMult;RefMult", 80, -0.5,
                           79.5, 80, -0.5, 79.5);

  hrefmult = new TH1F("hrefmult", ";RefMult;#it{N}_{events}", 80, -0.5, 79.5);
  hrefmult_vpdcut = new TH1F(
      "hrefmult_vpdcut", ";#it{N}_{refmult};#it{N}_{events}", 80, -0.5, 79.5);
  hrefmult_Ntofcut = new TH1F(
      "hrefmult_Ntofcut", ";#it{N}_{refmult};#it{N}_{events}", 80, -0.5, 79.5);
  hrefmult_Ntofcut_vpdcut =
      new TH1F("hrefmult_Ntof_vpdcut", ";#it{N}_{refmult};#it{N}_{events}", 80,
               -0.5, 79.5);

  hrefmult3 =
      new TH1F("hrefmult3", ";RefMult3;#it{N}_{events}", 80, -0.5, 79.5);
  hrefmult3_vpdcut = new TH1F(
      "hrefmult3_vpdcut", ";#it{N}_{refmult3};#it{N}_{events}", 80, -0.5, 79.5);
  hrefmult3_Ntofcut =
      new TH1F("hrefmult3_Ntofcut", ";#it{N}_{refmult3};#it{N}_{events}", 80,
               -0.5, 79.5);
  hrefmult3_Ntofcut_vpdcut =
      new TH1F("hrefmult3_Ntof_vpdcut", ";#it{N}_{refmult3};#it{N}_{events}",
               80, -0.5, 79.5);

#ifdef MINI_TREE
  mOutTree = new TTree("TrackInfo", "TrackInfo");
  mOutTree->Branch("Mom", mMom_Minitree, "Mom[3]/D");
  mOutTree->Branch("Charge", &mCharge_Minitree, "Charge/S");
  mOutTree->Branch("nSigmaProton", &mNSigmaProton_Minitree, "nSigmaProton/D");
  mOutTree->Branch("nSigmaKaon", &mNSigmaKaon_Minitree, "nSigmaKaon/D");
  mOutTree->Branch("nSigmaPion", &mNSigmaPion_Minitree, "nSigmaPion/D");
  mOutTree->Branch("BTofM2", &mBTofM2_Minitree_Minitree, "BTofM2/D");
#endif
  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Finish() {
  cout << "Writing histograms." << endl;

  mOutFile->cd();

  hnevents->Write();
  hvz->Write();
  hvzvzvpd->Write();
  hvzvzvpd_cut->Write();
  hbtofmatchmult->Write();
  h2NtofRefmult->Write();
  h2NtofRefmult_vpdcut->Write();
  hrefmult->Write();
  hrefmult_vpdcut->Write();
  hrefmult_Ntofcut->Write();
  hrefmult_Ntofcut_vpdcut->Write();
  hrefmult3->Write();
  hrefmult3_vpdcut->Write();
  hrefmult3_Ntofcut->Write();
  hrefmult3_Ntofcut_vpdcut->Write();
  mOutFile->Close();

#ifdef MINI_TREE
  mOutTree->Write();
  cout << mOutTree->GetEntries() << " tracks have been recorded." << endl;
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

  const TVector3 PrimaryVertex(mEvent->primaryVertex().x(),
                               mEvent->primaryVertex().y(),
                               mEvent->primaryVertex().z());

  hnevents->Fill(1.);
  hvz->Fill(PrimaryVertex.Z());
  hvzvzvpd->Fill(PrimaryVertex.Z(), mEvent->vzVpd());
  hbtofmatchmult->Fill(mEvent->nBTOFMatch());
  h2NtofRefmult->Fill(mEvent->nBTOFMatch(), mEvent->refMult());

  hrefmult->Fill(mEvent->refMult());
  hrefmult3->Fill(mEvent->refMult3());

  if (fabs(PrimaryVertex.Z() - mEvent->vzVpd()) < 6.) {
    hvzvzvpd_cut->Fill(PrimaryVertex.Z(), mEvent->vzVpd());
    h2NtofRefmult->Fill(mEvent->nBTOFMatch(), mEvent->refMult());
    hrefmult_vpdcut->Fill(mEvent->refMult());
    hrefmult3_vpdcut->Fill(mEvent->refMult3());
  }

  bool bool_ntofcut = true;

  if (bool_ntofcut) {
    hrefmult_Ntofcut->Fill(mEvent->refMult());
    hrefmult3_Ntofcut->Fill(mEvent->refMult3());
  }

  if (bool_ntofcut && fabs(PrimaryVertex.Z() - mEvent->vzVpd()) > 6.) {
    hrefmult_Ntofcut_vpdcut->Fill(mEvent->refMult());
    hrefmult3_Ntofcut_vpdcut->Fill(mEvent->refMult3());
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
  // if(TMath::Sqrt(TMath::Power(Vx-StMyCuts::VrCen[0], 2.)+TMath::Power(Vy-StMyCuts::VrCen[1],
  // 2.))>=StMyCuts::VrCut) 	return kFALSE;
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
