#include <map>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TString.h"
#include "TVector3.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StChain/StMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StThreeVector.hh"

#include "StMyCuts.h"
#include "StMyMaker.h"

ClassImp(StMyMaker)

//-----------------------------------------------------------------------------
StMyMaker::StMyMaker(const Char_t* const name, const Char_t* const inName, const Char_t* const outName, StPicoDstMaker* const picoDstMaker) : StMaker(name), mInName(inName), mOutName(outName), mPicoDstMaker(picoDstMaker){
	mOutFile = nullptr;
	Clear();
}

//-----------------------------------------------------------------------------
StMyMaker::~StMyMaker(){
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Init(){
	cout << "Initializing histograms." << endl;

	mOutFile = new TFile(mOutName, "recreate");
//	mOutFile->cd();

	const Int_t nRunIndices = StMyCuts::RunIdIndex.size();

	hnevents = new TH1D("hnevents", "", 2, -0.5, 1.5);
	hvxvy = new TH2F("hvxvy", ";#it{V}_{#it{x}} (cm);#it{V}_{#it{y}} (cm)", 100, -5., 5., 100, -5., 5.);
	hvz = new TH1F("hvz", ";#it{V}_{#it{z}} (cm);#it{N}_{events}", 80, 196., 204.);
	hvzvzvpd = new TH2F("hvzvzvpd", ";#it{V}_{#it{z}} (cm);#it{V}_{#it{z}}^{VPD} (cm)", 80, 196., 204., 80, 196., 204.);
	hnvpdhitseast = new TH1F("hnvpdhitseast", ";nVpdHitsEast;#it{N}_{events}", 50, -0.5, 49.5);
	hnvpdhitswest = new TH1F("hnvpdhitswest", ";nVpdHitsWest;#it{N}_{events}", 50, -0.5, 49.5);
	hbbcx = new TH1F("hbbcx", ";BBCx (Hz);#it{N}_{events}", 1000, -0.5, 99999.5);
	hzdcx = new TH1F("hzdcx", ";ZDCx (Hz);#it{N}_{events}", 1000, -0.5, 999.5);
	hrefmult = new TH1F("hrefmult", ";RefMult;#it{N}_{events}", 1000, -0.5, 999.5);
	hrefmult2 = new TH1F("hrefmult2", ";RefMult2;#it{N}_{events}", 1000, -0.5, 999.5);
	hrefmult3 = new TH1F("hrefmult3", ";RefMult3;#it{N}_{events}", 1000, -0.5, 999.5);
	hrefmult4 = new TH1F("hrefmult4", ";RefMult4;#it{N}_{events}", 1000, -0.5, 999.5);
	hgrefmult = new TH1F("hgrefmult", ";gRefMult;#it{N}_{events}", 1000, -0.5, 999.5);
	hfxtmult = new TH1F("hfxtmult", ";FxtMult;#it{N}_{events}", 1000, -0.5, 999.5);
	hfxtmult3 = new TH1F("hfxtmult3", ";FxtMult3;#it{N}_{events}", 1000, -0.5, 999.5);
	hfxtmult4 = new TH1F("hfxtmult4", ";FxtMult4;#it{N}_{events}", 1000, -0.5, 999.5);
	hbtofmatchmult = new TH1F("hbtofmatchmult", ";bTofMatchMult;#it{N}_{events}", 1000, -0.5, 999.5);
	hbtoftraymult = new TH1F("hbtoftraymult", ";bTofTrayMult;#it{N}_{events}", 1000, -0.5, 1999.5);
	hetofhitmult = new TH1F("hetofhitmult", ";eTofHitMult;#it{N}_{events}", 1000, -0.5, 999.5);
	hetofdigimult = new TH1F("hetofdigimult", ";eTofDigiMult;#it{N}_{events}", 1000, -0.5, 1999.5);
	hsumnmip = new TH1F("hsumnmip", ";#SigmanMIP;#it{N}_{events}", 1000, 0., 4000.);
	hsumtnmip = new TH1F("hsumtnmip", ";#Sigmatruncated nMIP;#it{N}_{events}", 1000, 0., 2000.);
	hnptracks = new TH1F("hnptracks", ";nPositiveTracks;#it{N}_{events}", 1000, -0.5, 999.5);
	hnntracks = new TH1F("hnntracks", ";nNegativeTracks;#it{N}_{events}", 1000, -0.5, 999.5);
	hnelectronp = new TH1F("hnelectronp", ";#it{n} e^{#plus};#it{N}_{events}", 500, -0.5, 499.5);
	hnelectronm = new TH1F("hnelectronm", ";#it{n} e^{#minus};#it{N}_{events}", 500, -0.5, 499.5);
	hnpionp = new TH1F("hnpionp", ";#it{n} #pi^{#plus};#it{N}_{events}", 1000, -0.5, 999.5);
	hnpionm = new TH1F("hnpionm", ";#it{n} #pi^{#minus};#it{N}_{events}", 1000, -0.5, 999.5);
	hnkaonp = new TH1F("hnkaonp", ";#it{n} K^{#plus};#it{N}_{events}", 500, -0.5, 499.5);
	hnkaonm = new TH1F("hnkaonm", ";#it{n} K^{#minus};#it{N}_{events}", 500, -0.5, 499.5);
	hnprotonp = new TH1F("hnprotonp", ";#it{n} p^{#plus};#it{N}_{events}", 500, -0.5, 499.5);
	hnprotonm = new TH1F("hnprotonm", ";#it{n} p^{#minus};#it{N}_{events}", 500, -0.5, 499.5);

	pvx = new TProfile("pvx", ";run index;<#it{V}_{#it{x}}>", nRunIndices, -0.5, nRunIndices-0.5);
	pvy = new TProfile("pvy", ";run index;<#it{V}_{#it{y}}>", nRunIndices, -0.5, nRunIndices-0.5);
	pvz = new TProfile("pvz", ";run index;<#it{V}_{#it{z}}>", nRunIndices, -0.5, nRunIndices-0.5);
	pvr = new TProfile("pvr", ";run index;<#it{V}_{#it{r}}>", nRunIndices, -0.5, nRunIndices-0.5);
	pvzmvzvpd = new TProfile("pvzmvzvpd", ";run index;<#it{V}_{#it{z}}#minus#it{V}_{#it{z}}^{VPD}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnvpdhitseast = new TProfile("pnvpdhitseast", ";run index;<nVpdHitsEast>", nRunIndices, -0.5, nRunIndices-0.5);
	pnvpdhitswest = new TProfile("pnvpdhitswest", ";run index;<nVpdHitsWest>", nRunIndices, -0.5, nRunIndices-0.5);
	pbbcx = new TProfile("pbbcx", ";run index;<BBCx>", nRunIndices, -0.5, nRunIndices-0.5);
	pzdcx = new TProfile("pzdcx", ";run index;<ZDCx>", nRunIndices, -0.5, nRunIndices-0.5);
	prefmult = new TProfile("prefmult", ";run index;<RefMult>", nRunIndices, -0.5, nRunIndices-0.5);
	prefmult2 = new TProfile("prefmult2", ";run index;<RefMult2>", nRunIndices, -0.5, nRunIndices-0.5);
	prefmult3 = new TProfile("prefmult3", ";run index;<RefMult3>", nRunIndices, -0.5, nRunIndices-0.5);
	prefmult4 = new TProfile("prefmult4", ";run index;<RefMult4>", nRunIndices, -0.5, nRunIndices-0.5);
	pgrefmult = new TProfile("pgrefmult", ";run index;<gRefMult>", nRunIndices, -0.5, nRunIndices-0.5);
	pfxtmult = new TProfile("pfxtmult", ";run index;<FxtMult>", nRunIndices, -0.5, nRunIndices-0.5);
	pfxtmult3 = new TProfile("pfxtmult3", ";run index;<FxtMult3>", nRunIndices, -0.5, nRunIndices-0.5);
	pfxtmult4 = new TProfile("pfxtmult4", ";run index;<FxtMult4>", nRunIndices, -0.5, nRunIndices-0.5);
	pbtofmatchmult = new TProfile("pbtofmatchmult", ";run index;<bTofMatchMult>", nRunIndices, -0.5, nRunIndices-0.5);
	pbtoftraymult = new TProfile("pbtoftraymult", ";run index;<bTofTrayMult>", nRunIndices, -0.5, nRunIndices-0.5);
	petofhitmult = new TProfile("petofhitmult", ";run index;<eTofHitMult>", nRunIndices, -0.5, nRunIndices-0.5);
	petofdigimult = new TProfile("petofdigimult", ";run index;<eTofDigiMult>", nRunIndices, -0.5, nRunIndices-0.5);
	psumnmip = new TProfile("psumnmip", ";run index;<#SigmanMIP>", nRunIndices, -0.5, nRunIndices-0.5);
	psumtnmip = new TProfile("psumtnmip", ";run index;<#Sigmatruncated nMIP>", nRunIndices, -0.5, nRunIndices-0.5);
	pnptracks = new TProfile("pnptracks", ";run index;<nPositiveTracks>", nRunIndices, -0.5, nRunIndices-0.5);
	pnntracks = new TProfile("pnntracks", ";run index;<nNegativeTracks>", nRunIndices, -0.5, nRunIndices-0.5);
	pnelectronp = new TProfile("pnelectronp", ";run index;<#it{n} e^{#plus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnelectronm = new TProfile("pnelectronm", ";run index;<#it{n} e^{#minus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnpionp = new TProfile("pnpionp", ";run index;<#it{n} #pi^{#plus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnpionm = new TProfile("pnpionm", ";run index;<#it{n} #pi^{#minus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnkaonp = new TProfile("pnkaonp", ";run index;<#it{n} K^{#plus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnkaonm = new TProfile("pnkaonm", ";run index;<#it{n} K^{#minus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnprotonp = new TProfile("pnprotonp", ";run index;<#it{n} p^{#plus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnprotonm = new TProfile("pnprotonm", ";run index;<#it{n} p^{#minus}>", nRunIndices, -0.5, nRunIndices-0.5);
	pq1xtpc = new TProfile("pq1xtpc", ";run index;<#it{Q}_{1#it{x}}^{TPC}>", nRunIndices, -0.5, nRunIndices-0.5);
	pq1ytpc = new TProfile("pq1ytpc", ";run index;<#it{Q}_{1#it{y}}^{TPC}>", nRunIndices, -0.5, nRunIndices-0.5);
	pq2xtpc = new TProfile("pq2xtpc", ";run index;<#it{Q}_{2#it{x}}^{TPC}>", nRunIndices, -0.5, nRunIndices-0.5);
	pq2ytpc = new TProfile("pq2ytpc", ";run index;<#it{Q}_{2#it{y}}^{TPC}>", nRunIndices, -0.5, nRunIndices-0.5);
	
	hdca = new TH1F("hdca", ";Dca (cm);#it{N}_{tracks}", 400, 0., 4.);
	hdcaxy = new TH1F("hdcaxy", ";Dca_{#it{xy}} (cm);#it{N}_{tracks}", 400, 0., 4.);
	hdcaphi = new TH1F("hdcaphi", ";Dca_{#it{#varphi}};#it{N}_{tracks}", 360, -TMath::Pi(), TMath::Pi());
	hdcaz = new TH1F("hdcaz", ";Dca_{#it{z}} (cm);#it{N}_{tracks}", 400, 0., 4.);
	hsdcaxy = new TH1F("hsdcaxy", ";signed Dca_{#it{xy}} (cm);#it{N}_{tracks}", 800, -4., 4.);
	hnhitsfit = new TH1F("hnhitsfit", ";nHitsFit;#it{N}_{tracks}", 100, -0.5, 99.5);
	hnhitsratio = new TH1F("hnhitsratio", ";nHitsFit/nHitsPoss;#it{N}_{tracks}", 110, 0., 1.1);
	hnhitsdedx = new TH1F("hnhitsdedx", ";nHitsDedx;#it{N}_{tracks}", 100, -0.5, 99.5);
	hpt = new TH1F("hpt", ";#it{p}_{T} (GeV/#it{c});#it{N}_{tracks}", 500, 0., 10.);
	hphi = new TH1F("hphi", ";#it{#varphi};#it{N}_{tracks}", 360, -TMath::Pi(), TMath::Pi());
	heta = new TH1F("heta", ";#it{#eta};#it{N}_{tracks}", 300, -3., 0.);
	hpqdedx = new TH2F("hpqdedx", ";#it{p}#it{q} (GeV/#it{c});d#it{E}/d#it{x} (keV/cm)", 200, -10., 10., 200, 0., 20.);
	hpqnsigmaelectron = new TH2F("hpqnsigmaelectron", ";#it{p}#it{q} (GeV/#it{c});#it{n}#it{#sigma}_{electron}", 200, -10., 10., 200, -20., 20.);
	hpqnsigmapion = new TH2F("hpqnsigmapion", ";#it{p}#it{q} (GeV/#it{c});#it{n}#it{#sigma}_{pion}", 200, -10., 10., 200, -20., 20.);
	hpqnsigmakaon = new TH2F("hpqnsigmakaon", ";#it{p}#it{q} (GeV/#it{c});#it{n}#it{#sigma}_{kaon}", 200, -10., 10., 200, -20., 20.);
	hpqnsigmaproton = new TH2F("hpqnsigmaproton", ";#it{p}#it{q} (GeV/#it{c});#it{n}#it{#sigma}_{proton}", 200, -10., 10., 200, -20., 20.);
	hbtofyzlocal = new TH2F("hbtofyzlocal", ";bTofYLocal (cm);bTofZLocal (cm)", 120, -6., 6., 80, -4., 4.);
	hbtof1obetacheck = new TH2F("hbtof1obetacheck", ";bTOF 1/#it{#beta};bTOF 1/#it{#beta} (recalculation)", 200, 0.7, 2.7, 200, 0.7, 2.7);
	hpqbtof1obeta = new TH2F("hpqbtof1obeta", ";#it{p}#it{q} (GeV/#it{c});bTOF 1/#it{#beta}", 200, -10., 10., 200, 0.7, 2.7);
	hpqbtof1obetare = new TH2F("hpqbtof1obetare", ";#it{p}#it{q} (GeV/#it{c});bTOF 1/#it{#beta} (recalculation)", 200, -10., 10., 200, 0.7, 2.7);
	hpqbtofm2 = new TH2F("hpqbtofm2", ";#it{p}#it{q} (GeV/#it{c});bTOF #it{m}^{2} (GeV^{2}/#it{c}^{4})", 200, -10., 10., 250, -0.2, 2.3);
	hpqbtofm2re = new TH2F("hpqbtofm2re", ";#it{p}#it{q} (GeV/#it{c});bTOF #it{m}^{2} (recalculation) (GeV^{2}/#it{c}^{4})", 200, -10., 10., 250, -0.2, 2.3);
	hetofdeltaxy = new TH2F("hetofdeltaxy", ";eTofDeltaX (cm);eTofDeltaY (cm)", 120, -6., 6., 240, -12., 12.);
	hetof1obetacheck = new TH2F("hetof1obetacheck", ";eTOF 1/#it{#beta};eTOF 1/#it{#beta} (recalculation)", 200, 0.7, 2.7, 200, 0.7, 2.7);
	hpqetof1obeta = new TH2F("hpqetof1obeta", ";#it{p}#it{q} (GeV/#it{c});eTOF 1/#it{#beta}", 200, -10., 10., 200, 0.7, 2.7);
	hpqetof1obetare = new TH2F("hpqetof1obetare", ";#it{p}#it{q} (GeV/#it{c});eTOF 1/#it{#beta} (recalculation)", 200, -10., 10., 200, 0.7, 2.7);
	hpqetofm2 = new TH2F("hpqetofm2", ";#it{p}#it{q} (GeV/#it{c});eTOF #it{m}^{2} (GeV^{2}/#it{c}^{4})", 200, -10., 10., 250, -0.2, 2.3);
	hpqetofm2re = new TH2F("hpqetofm2re", ";#it{p}#it{q} (GeV/#it{c});eTOF #it{m}^{2} (recalculation) (GeV^{2}/#it{c}^{4})", 200, -10., 10., 250, -0.2, 2.3);

	pdca = new TProfile("pdca", ";run index;<Dca>", nRunIndices, -0.5, nRunIndices-0.5);
	pdcaxy = new TProfile("pdcaxy", ";run index;<Dca_{#it{xy}}>", nRunIndices, -0.5, nRunIndices-0.5);
	pdcaphi = new TProfile("pdcaphi", ";run index;<Dca #it{#varphi}>", nRunIndices, -0.5, nRunIndices-0.5);
	pdcaz = new TProfile("pdcaz", ";run index;<Dca_{#it{z}}>", nRunIndices, -0.5, nRunIndices-0.5);
	psdcaxy = new TProfile("psdcaxy", ";run index;<signed Dca_{#it{xy}}>", nRunIndices, -0.5, nRunIndices-0.5);
	psdcaxy2 = new TProfile("psdcaxy2", ";run index;<(signed Dca_{#it{xy}})^{2}>", nRunIndices, -0.5, nRunIndices-0.5);
	pnhitsfit = new TProfile("pnhitsfit", ";run index;<nHitsFit>", nRunIndices, -0.5, nRunIndices-0.5);
	pnhitsratio = new TProfile("pnhitsratio", ";run index;<nHitsFit/nHitsPoss>", nRunIndices, -0.5, nRunIndices-0.5);
	pnhitsdedx = new TProfile("pnhitsdedx", ";run index;<nHitsDedx>", nRunIndices, -0.5, nRunIndices-0.5);
	ppt = new TProfile("ppt", ";run index;<#it{p}_{T}>", nRunIndices, -0.5, nRunIndices-0.5);
	pphi = new TProfile("pphi", ";run index;<#it{#varphi}>", nRunIndices, -0.5, nRunIndices-0.5);
	peta = new TProfile("peta", ";run index;<#it{#eta}>", nRunIndices, -0.5, nRunIndices-0.5);
	pdedx = new TProfile("pdedx", ";run index;<d#it{E}/d#it{x}>", nRunIndices, -0.5, nRunIndices-0.5);
	pbtof1obeta = new TProfile("pbtof1obeta", ";run index;<bTOF 1/#it{#beta}>", nRunIndices, -0.5, nRunIndices-0.5);
	pbtofylocal = new TProfile("pbtofylocal", ";run index;<bTofYLocal>", nRunIndices, -0.5, nRunIndices-0.5);
	pbtofzlocal = new TProfile("pbtofzlocal", ";run index;<bTofZLocal>", nRunIndices, -0.5, nRunIndices-0.5);
	petof1obeta = new TProfile("petof1obeta", ";run index;<eTOF 1/#it{#beta}>", nRunIndices, -0.5, nRunIndices-0.5);
	petofdeltax = new TProfile("petofdeltax", ";run index;<eTofDeltaX>", nRunIndices, -0.5, nRunIndices-0.5);
	petofdeltay = new TProfile("petofdeltay", ";run index;<eTofDeltaY>", nRunIndices, -0.5, nRunIndices-0.5);
	
	return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Finish(){
	cout << "Writing histograms." << endl;

	mOutFile->cd();

	hnevents->Write();
	hvxvy->Write();
	hvz->Write();
	hvzvzvpd->Write();
	hnvpdhitseast->Write();
	hnvpdhitswest->Write();
	hbbcx->Write();
	hzdcx->Write();
	hrefmult->Write();
	hrefmult2->Write();
	hrefmult3->Write();
	hrefmult4->Write();
	hgrefmult->Write();
	hfxtmult->Write();
	hfxtmult3->Write();
	hfxtmult4->Write();
	hbtofmatchmult->Write();
	hbtoftraymult->Write();
	hetofhitmult->Write();
	hetofdigimult->Write();
	hsumnmip->Write();
	hsumtnmip->Write();
	hnptracks->Write();
	hnntracks->Write();
	hnelectronp->Write();
	hnelectronm->Write();
	hnpionp->Write();
	hnpionm->Write();
	hnkaonp->Write();
	hnkaonm->Write();
	hnprotonp->Write();
	hnprotonm->Write();

	pvx->Write();
	pvy->Write();
	pvz->Write();
	pvr->Write();
	pvzmvzvpd->Write();
	pnvpdhitseast->Write();
	pnvpdhitswest->Write();
	pbbcx->Write();
	pzdcx->Write();
	prefmult->Write();
	prefmult2->Write();
	prefmult3->Write();
	prefmult4->Write();
	pgrefmult->Write();
	pfxtmult->Write();
	pfxtmult3->Write();
	pfxtmult4->Write();
	pbtofmatchmult->Write();
	pbtoftraymult->Write();
	petofhitmult->Write();
	petofdigimult->Write();
	psumnmip->Write();
	psumtnmip->Write();
	pnptracks->Write();
	pnntracks->Write();
	pnelectronp->Write();
	pnelectronm->Write();
	pnpionp->Write();
	pnpionm->Write();
	pnkaonp->Write();
	pnkaonm->Write();
	pnprotonp->Write();
	pnprotonm->Write();
	pq1xtpc->Write();
	pq1ytpc->Write();
	pq2xtpc->Write();
	pq2ytpc->Write();

	hdca->Write();
	hdcaxy->Write();
	hdcaphi->Write();
	hdcaz->Write();
	hsdcaxy->Write();
	hnhitsfit->Write();
	hnhitsratio->Write();
	hnhitsdedx->Write();
	hpt->Write();
	hphi->Write();
	heta->Write();
	hpqdedx->Write();
	hpqnsigmaelectron->Write();
	hpqnsigmapion->Write();
	hpqnsigmakaon->Write();
	hpqnsigmaproton->Write();
	hbtofyzlocal->Write();
	hbtof1obetacheck->Write();
	hpqbtof1obeta->Write();
	hpqbtof1obetare->Write();
	hpqbtofm2->Write();
	hpqbtofm2re->Write();
	hetofdeltaxy->Write();
	hetof1obetacheck->Write();
	hpqetof1obeta->Write();
	hpqetof1obetare->Write();
	hpqetofm2->Write();
	hpqetofm2re->Write();

	pdca->Write();
	pdcaxy->Write();
	pdcaphi->Write();
	pdcaz->Write();
	psdcaxy->Write();
	psdcaxy2->Write();
	pnhitsfit->Write();
	pnhitsratio->Write();
	pnhitsdedx->Write();
	ppt->Write();
	pphi->Write();
	peta->Write();
	pdedx->Write();
	pbtof1obeta->Write();
	pbtofylocal->Write();
	pbtofzlocal->Write();
	petof1obeta->Write();
	petofdeltax->Write();
	petofdeltay->Write();

	mOutFile->Close();

	return kStOK;
}

//-----------------------------------------------------------------------------
void StMyMaker::Clear(Option_t* option){
	mPicoDst       = nullptr;
	mEvent         = nullptr;
	mTrack         = nullptr;
	mBTofPidTraits = nullptr;
	mETofPidTraits = nullptr;
}

//-----------------------------------------------------------------------------
Int_t StMyMaker::Make(){
	if(!mPicoDstMaker){
		LOG_WARN << "No PicoDstMaker! Skip!" << endm;
		return kStWarn;
	}

	mPicoDst = (const StPicoDst*)mPicoDstMaker->picoDst();
	if(!mPicoDst){
		LOG_WARN << "No PicoDst! Skip!" << endm;
		return kStWarn;
	}
	
	MakeEvent();

	return kStOK;
}

//-----------------------------------------------------------------------------
const Int_t StMyMaker::MakeEvent(){
	mEvent = (const StPicoEvent*)mPicoDst->event();
	if(!mEvent){
		LOG_WARN << "Error opening picoDst Event! Skip!" << endm;
		return kStWarn;
	}

	hnevents->Fill(0.);

	if(!isGoodTrigger() || !isGoodEvent())
		return kStOK;

	mRunIndex = StMyCuts::RunIdIndex.at(mEvent->runId());
	const TVector3 PrimaryVertex = mEvent->primaryVertex();

	hnevents->Fill(1.);
	hvxvy->Fill(PrimaryVertex.X(), PrimaryVertex.Y());
	hvz->Fill(PrimaryVertex.Z());
	hvzvzvpd->Fill(PrimaryVertex.Z(), mEvent->vzVpd());
	hnvpdhitseast->Fill(mEvent->nVpdHitsEast());
	hnvpdhitswest->Fill(mEvent->nVpdHitsWest());
	hbbcx->Fill(mEvent->BBCx());
	hzdcx->Fill(mEvent->ZDCx());
	hrefmult->Fill(mEvent->refMult());
	hrefmult2->Fill(mEvent->refMult2());
	hrefmult3->Fill(mEvent->refMult3());
	hrefmult4->Fill(mEvent->refMult4());
	hgrefmult->Fill(mEvent->grefMult());
	hfxtmult->Fill(mEvent->fxtMult());
	hbtofmatchmult->Fill(mEvent->nBTOFMatch());
	hbtoftraymult->Fill(mEvent->btofTrayMultiplicity());
	hetofhitmult->Fill(mEvent->etofHitMultiplicity());
	hetofdigimult->Fill(mEvent->etofDigiMultiplicity());

	pvx->Fill(mRunIndex, PrimaryVertex.X());
	pvy->Fill(mRunIndex, PrimaryVertex.Y());
	pvz->Fill(mRunIndex, PrimaryVertex.Z());
	pvr->Fill(mRunIndex, TMath::Sqrt(TMath::Power(PrimaryVertex.X()-StMyCuts::VrCen[0], 2.)+TMath::Power(PrimaryVertex.Y()-StMyCuts::VrCen[1], 2.)));
	pvzmvzvpd->Fill(mRunIndex, PrimaryVertex.Z()-mEvent->vzVpd());
	pnvpdhitseast->Fill(mRunIndex, mEvent->nVpdHitsEast());
	pnvpdhitswest->Fill(mRunIndex, mEvent->nVpdHitsWest());
	pbbcx->Fill(mRunIndex, mEvent->BBCx());
	pzdcx->Fill(mRunIndex, mEvent->ZDCx());
	prefmult->Fill(mRunIndex, mEvent->refMult());
	prefmult2->Fill(mRunIndex, mEvent->refMult2());
	prefmult3->Fill(mRunIndex, mEvent->refMult3());
	prefmult4->Fill(mRunIndex, mEvent->refMult4());
	pgrefmult->Fill(mRunIndex, mEvent->grefMult());
	pfxtmult->Fill(mRunIndex, mEvent->fxtMult());
	pbtofmatchmult->Fill(mRunIndex, mEvent->nBTOFMatch());
	pbtoftraymult->Fill(mRunIndex, mEvent->btofTrayMultiplicity());
	petofhitmult->Fill(mRunIndex, mEvent->etofHitMultiplicity());
	petofdigimult->Fill(mRunIndex, mEvent->etofDigiMultiplicity());

	mFxtMult3 = mFxtMult4 = mNPTracks = mNNTracks = mNElectronP = mNElectronM = mNPionP = mNPionM = mNKaonP = mNKaonM = mNProtonP = mNProtonM = 0;
	mQ1xTpc = mQ1yTpc = mQ2xTpc = mQ2yTpc = 0.;
	const Int_t nTracks = mPicoDst->numberOfTracks();
	for(Int_t it=0; it<nTracks; it++)
		MakeTrack(it);
	
	hfxtmult3->Fill(mFxtMult3);
	hfxtmult4->Fill(mFxtMult4);
	hnptracks->Fill(mNPTracks);
	hnntracks->Fill(mNNTracks);
	hnelectronp->Fill(mNElectronP);
	hnelectronm->Fill(mNElectronM);
	hnpionp->Fill(mNPionP);
	hnpionm->Fill(mNPionM);
	hnkaonp->Fill(mNKaonP);
	hnkaonm->Fill(mNKaonM);
	hnprotonp->Fill(mNProtonP);
	hnprotonm->Fill(mNProtonM);

	pfxtmult3->Fill(mRunIndex, mFxtMult3);
	pfxtmult4->Fill(mRunIndex, mFxtMult4);
	pnptracks->Fill(mRunIndex, mNPTracks);
	pnntracks->Fill(mRunIndex, mNNTracks);
	pnelectronp->Fill(mRunIndex, mNElectronP);
	pnelectronm->Fill(mRunIndex, mNElectronM);
	pnpionp->Fill(mRunIndex, mNPionP);
	pnpionm->Fill(mRunIndex, mNPionM);
	pnkaonp->Fill(mRunIndex, mNKaonP);
	pnkaonm->Fill(mRunIndex, mNKaonM);
	pnprotonp->Fill(mRunIndex, mNProtonP);
	pnprotonm->Fill(mRunIndex, mNProtonM);
	pq1xtpc->Fill(mRunIndex, mQ1xTpc);
	pq1ytpc->Fill(mRunIndex, mQ1yTpc);
	pq2xtpc->Fill(mRunIndex, mQ2xTpc);
	pq2ytpc->Fill(mRunIndex, mQ2yTpc);

	hsumnmip->Fill(mSumNMip);
	hsumtnmip->Fill(mSumTNMip);

	psumnmip->Fill(mRunIndex, mSumNMip);
	psumtnmip->Fill(mRunIndex, mSumTNMip);
	
	return kStOK;
}

//-----------------------------------------------------------------------------
const Int_t StMyMaker::MakeTrack(const Int_t it){
	mTrack = (const StPicoTrack*)mPicoDst->track(it);
	if(!mTrack){
		LOG_WARN << "Error opening picoDst Track! Skip!" << endm;
		return kStWarn;
	}
	mBTofPidTraits = mTrack->isTofTrack()?(const StPicoBTofPidTraits*)mPicoDst->btofPidTraits(mTrack->bTofPidTraitsIndex()):nullptr;
	mETofPidTraits = mTrack->isETofTrack()?(const StPicoETofPidTraits*)mPicoDst->etofPidTraits(mTrack->eTofPidTraitsIndex()):nullptr;

	if(mTrack->isPrimary() && mTrack->nHitsFit()>10 && mTrack->nSigmaProton()<-3.)
		mFxtMult3 ++;
	if(mTrack->isPrimary() && mTrack->nHitsFit()>10 && TMath::Abs(mTrack->nSigmaKaon())>3.)
		mFxtMult4 ++;

	if(!isGoodTrack())
		return kStOK;
	
	const TVector3 Dca = mTrack->gDCA(mEvent->primaryVertex());
	const Double_t SDcaXY = mTrack->gDCAs(mEvent->primaryVertex());
	const TVector3 PMom = mTrack->pMom();
	const Double_t PPt = PMom.Pt();
	const Double_t PPhi= PMom.Phi();
	const Double_t PEta = PMom.Eta();
	const Double_t PP2 = PMom.Mag2();
	const Double_t PPQ = PMom.Mag()*mTrack->charge();
	const Double_t BTofBeta = getBTofBeta(0);
	const Double_t BTofBetaRe = getBTofBeta(2);
	const Double_t BTofM2 = BTofBeta<1.e-5?-999.:PP2*(TMath::Power(BTofBeta, -2.)-1.);
	const Double_t BTofM2Re = BTofBetaRe<1.e-5?-999.:PP2*(TMath::Power(BTofBetaRe, -2.)-1.);
	const Double_t ETofBeta = getETofBeta(0);
	const Double_t ETofBetaRe = getETofBeta(2);
	const Double_t ETofM2 = ETofBeta<1.e-5?-999.:PP2*(TMath::Power(ETofBeta, -2.)-1.);
	const Double_t ETofM2Re = ETofBetaRe<1.e-5?-999.:PP2*(TMath::Power(ETofBetaRe, -2.)-1.);
	
	hdca->Fill(Dca.Mag());
	hdcaxy->Fill(Dca.Pt());
	hdcaphi->Fill(Dca.Phi());
	hdcaz->Fill(Dca.Z());
	hsdcaxy->Fill(mTrack->gDCAs(mEvent->primaryVertex()));
	hnhitsfit->Fill(mTrack->nHitsFit());
	hnhitsratio->Fill((const Double_t)mTrack->nHitsFit()/mTrack->nHitsPoss());
	hnhitsdedx->Fill(mTrack->nHitsDedx());
	hpt->Fill(PPt);
	hphi->Fill(PPhi);
	heta->Fill(PEta);
	hpqdedx->Fill(PPQ, mTrack->dEdx());
	hpqnsigmaelectron->Fill(PPQ, mTrack->nSigmaElectron());
	hpqnsigmapion->Fill(PPQ, mTrack->nSigmaPion());
	hpqnsigmakaon->Fill(PPQ, mTrack->nSigmaKaon());
	hpqnsigmaproton->Fill(PPQ, mTrack->nSigmaProton());
	if(mTrack->isTofTrack() && mBTofPidTraits && mBTofPidTraits->btofMatchFlag()>0){
		hbtofyzlocal->Fill(mBTofPidTraits->btofYLocal(), mBTofPidTraits->btofZLocal());
		hbtof1obetacheck->Fill(1./BTofBeta, 1./BTofBetaRe);
	}
	if(!(BTofBeta<1.e-5)){
		hpqbtof1obeta->Fill(PPQ, 1./BTofBeta);
		hpqbtofm2->Fill(PPQ, BTofM2);
	}
	if(!(BTofBetaRe<1.e-5)){
		hpqbtof1obetare->Fill(PPQ, 1./BTofBetaRe);
		hpqbtofm2re->Fill(PPQ, BTofM2Re);
	}
	if(mTrack->isETofTrack() && mETofPidTraits && mETofPidTraits->matchFlag()>0){
		hetofdeltaxy->Fill(mETofPidTraits->deltaX(), mETofPidTraits->deltaY());
		hetof1obetacheck->Fill(1./ETofBeta, 1./ETofBetaRe);
	}
	if(!(ETofBeta<1.e-5)){
		hpqetof1obeta->Fill(PPQ, 1./ETofBeta);
		hpqetofm2->Fill(PPQ, ETofM2);
	}
	if(!(ETofBetaRe<1.e-5)){
		hpqetof1obetare->Fill(PPQ, 1./ETofBetaRe);
		hpqetofm2re->Fill(PPQ, ETofM2Re);
	}

	pdca->Fill(mRunIndex, Dca.Mag());
	pdcaxy->Fill(mRunIndex, Dca.Pt());
	pdcaphi->Fill(mRunIndex, Dca.Phi());
	pdcaz->Fill(mRunIndex, Dca.Z());
	psdcaxy->Fill(mRunIndex, SDcaXY);
	psdcaxy2->Fill(mRunIndex, TMath::Power(SDcaXY, 2.));
	pnhitsfit->Fill(mRunIndex, mTrack->nHitsFit());
	pnhitsratio->Fill(mRunIndex, (const Double_t)mTrack->nHitsFit()/mTrack->nHitsPoss());
	pnhitsdedx->Fill(mRunIndex, mTrack->nHitsDedx());
	ppt->Fill(mRunIndex, PPt);
	pphi->Fill(mRunIndex, PPhi);
	peta->Fill(mRunIndex, PEta);
	pdedx->Fill(mRunIndex, mTrack->dEdx());
	if(!(BTofBeta<1.e-5)){
		pbtof1obeta->Fill(mRunIndex, 1./BTofBeta);
		pbtofylocal->Fill(mRunIndex, mBTofPidTraits->btofYLocal());
		pbtofzlocal->Fill(mRunIndex, mBTofPidTraits->btofZLocal());
	}
	if(!(ETofBeta<1.e-5)){
		petof1obeta->Fill(mRunIndex, 1./ETofBeta);
		petofdeltax->Fill(mRunIndex, mETofPidTraits->deltaX());
		petofdeltay->Fill(mRunIndex, mETofPidTraits->deltaY());
	}

	if(mTrack->charge()>0){
		mNPTracks ++;
		if(isElectron())
			mNElectronP ++;
		if(isPion())
			mNPionP ++;
		if(isKaon())
			mNKaonP ++;
		if(isProton())
			mNProtonP ++;
	}
	if(mTrack->charge()<0){
		mNNTracks ++;
		if(isElectron())
			mNElectronM ++;
		if(isPion())
			mNPionM ++;
		if(isKaon())
			mNKaonM ++;
		if(isProton())
			mNProtonM ++;
	}

	if((PPt>StMyCuts::PtCutQ[0] && PPt<StMyCuts::PtCutQ[1]) && ((PEta>StMyCuts::EtaCutQ[0] && PEta<StMyCuts::EtaCutQ[1]) || (PEta>StMyCuts::EtaCutQ[2] && PEta<StMyCuts::EtaCutQ[3]))){
		const Double_t Weight = PPt<0.2?0.:TMath::Min(PPt, 2.);
		mQ1xTpc += -PEta*Weight*TMath::Cos(PPhi);
		mQ1yTpc += -PEta*Weight*TMath::Sin(PPhi);
		mQ2xTpc += Weight*TMath::Cos(2.*PPhi);
		mQ2yTpc += Weight*TMath::Sin(2.*PPhi);
	}

	return kStOK;
}

//-----------------------------------------------------------------------------
const Int_t StMyMaker::MakeEpdHit(const Int_t it){
	mEpdHit = (const StPicoEpdHit*)mPicoDst->epdHit(it);
	if(!mEpdHit){
		LOG_WARN << "Error opening picoDst EpdHit! Skip!" << endm;
		return kStWarn;
	}

	if(mEpdGeom->IsEast(mEpdHit->id()))
		mNEpdHitsEast ++;
	if(mEpdGeom->IsWest(mEpdHit->id()))
		mNEpdHitsWest ++;
	
	if(!mEpdHit->isGood())
		return kStOK;

	mSumNMip += mEpdHit->nMIP();
	mSumTNMip += mEpdHit->TnMIP(2., 0.3);

	const TVector3 StraightLine = mEpdGeom->TileCenter(mEpdHit->id())-mEvent->primaryVertex();
	const Double_t Weight = mEpdHit->TnMIP(2., 0.3);
	mQ1xEpd += StraightLine.Eta()*Weight*TMath::Cos(StraightLine.Phi());
	mQ1yEpd += StraightLine.Eta()*Weight*TMath::Sin(StraightLine.Phi());
	mQ2xEpd += Weight*TMath::Cos(2.*StraightLine.Phi());
	mQ2yEpd += Weight*TMath::Sin(2.*StraightLine.Phi());

	return kStOK;
}

//-----------------------------------------------------------------------------
const Double_t StMyMaker::getBTofBeta(const UChar_t iCase) const{
	if(!mTrack->isTofTrack() || !mBTofPidTraits || mBTofPidTraits->btofMatchFlag()<=0)
		return 0.;
	Double_t BTofBeta = mBTofPidTraits->btofBeta();
	if(iCase==0)
		return BTofBeta<1.e-5?0.:BTofBeta;
	if(iCase==1 && !(BTofBeta<1.e-5))
		return BTofBeta;
	const TVector3 Origin = mTrack->origin();
	const TVector3 BTofHitPos = mBTofPidTraits->btofHitPos();
	BTofBeta = tofPathLength(new StThreeVector<Double_t>(Origin.X(), Origin.Y(), Origin.Z()), new StThreeVector<Double_t>(BTofHitPos.X(), BTofHitPos.Y(), BTofHitPos.Z()), mTrack->helix(mEvent->bField()).curvature())/(mBTofPidTraits->btof()*C_C_LIGHT/1.e9);
	return BTofBeta<1.e-5?0.:BTofBeta;
}

//-----------------------------------------------------------------------------
const Double_t StMyMaker::getETofBeta(const UChar_t iCase) const{
	if(!mTrack->isETofTrack() || !mETofPidTraits || mETofPidTraits->matchFlag()<=0)
		return 0.;
	Double_t ETofBeta = mETofPidTraits->beta();
	if(iCase==0)
		return ETofBeta<1.e-5?0.:ETofBeta;
	if(iCase==1 && !(ETofBeta<1.e-5))
		return ETofBeta;
	const TVector3 Origin = mTrack->origin();
	const TVector3 ETofCrossingPos = mETofPidTraits->crossingPos();
	ETofBeta = tofPathLength(new StThreeVector<Double_t>(Origin.X(), Origin.Y(), Origin.Z()), new StThreeVector<Double_t>(ETofCrossingPos.X(), ETofCrossingPos.Y(), ETofCrossingPos.Z()), mTrack->helix(mEvent->bField()).curvature())/(mETofPidTraits->tof()*C_C_LIGHT/1.e9);
	return ETofBeta<1.e-5?0.:ETofBeta;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodTrigger() const{
	for(auto iTrigger : StMyCuts::Trigger)
		if(mEvent->isTrigger(iTrigger))
			return kTRUE;
	return kFALSE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodEvent() const{
	const Double_t Vx = mEvent->primaryVertex().X();
	const Double_t Vy = mEvent->primaryVertex().Y();
	const Double_t Vz = mEvent->primaryVertex().Z();
	if(TMath::Abs(Vx)<1.e-5 && TMath::Abs(Vy)<1.e-5 && TMath::Abs(Vz)<1.e-5)
		return kFALSE;
	if(Vz<=StMyCuts::VzCut[0] || Vz>=StMyCuts::VzCut[1])
		return kFALSE;
	if(TMath::Sqrt(TMath::Power(Vx-StMyCuts::VrCen[0], 2.)+TMath::Power(Vy-StMyCuts::VrCen[1], 2.))>=StMyCuts::VrCut)
		return kFALSE;
	return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isGoodTrack() const{
	if(!mTrack->isPrimary())
		return kFALSE;
	if(mTrack->gDCA(mEvent->primaryVertex()).Mag()>=StMyCuts::DcaCut)
		return kFALSE;
	if(mTrack->nHitsFit()<=StMyCuts::NHitsFitCut)
		return kFALSE;
	if((const Double_t)mTrack->nHitsFit()/mTrack->nHitsPoss()<=StMyCuts::NHitsRatioCut)
		return kFALSE;
	if(mTrack->nHitsDedx()<=StMyCuts::NHitsDedxCut)
		return kFALSE;
	return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isElectron() const{
	if(TMath::Abs(mTrack->nSigmaElectron())>=StMyCuts::NSigmaCut)
		return kFALSE;
	return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isPion() const{
	if(TMath::Abs(mTrack->nSigmaPion())>=StMyCuts::NSigmaCut)
		return kFALSE;
	return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isKaon() const{
	if(TMath::Abs(mTrack->nSigmaKaon())>=StMyCuts::NSigmaCut)
		return kFALSE;
	return kTRUE;
}

//-----------------------------------------------------------------------------
const Bool_t StMyMaker::isProton() const{
	if(TMath::Abs(mTrack->nSigmaProton())>=StMyCuts::NSigmaCut)
		return kFALSE;
	return kTRUE;
}