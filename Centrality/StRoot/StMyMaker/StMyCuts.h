#ifndef StMyCuts_h
#define StMyCuts_h

#include <map>
#include <vector>

#include "Rtypes.h"

namespace StMyCuts {
const std::vector<UInt_t> Trigger = {500904};
const std::vector<Double_t> VzCut = {-30., 30.};
// const Double_t VrCut = 2.;
const Double_t DcaCut = 1.;
const Int_t NHitsFitCut = 20;
const Double_t NHitsRatioCut = 0.52;
// const Int_t NHitsDedxCut = 5;
const Double_t NSigmaCut = 2.;
const std::vector<Double_t> PtCutQ = {0.2, 2.};
const std::vector<Double_t> EtaCutQ = {-1., 1.};

} // namespace StMyCuts

#endif