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

const std::map<Int_t, Int_t> RunIdIndex = {
    {16124017, 0},  {16124018, 1},  {16124019, 2},  {16124033, 3},
    {16124034, 4},  {16125003, 5},  {16125014, 6},  {16125015, 7},
    {16125016, 8},  {16125024, 9},  {16125035, 10}, {16125038, 11},
    {16125039, 12}, {16125042, 13}, {16125043, 14}, {16125046, 15},
    {16125050, 16}, {16125052, 17}, {16125053, 18}, {16125054, 19},
    {16125056, 20}, {16125057, 21}, {16125058, 22}, {16125060, 23},
    {16125061, 24}, {16125062, 25}, {16125063, 26}, {16126006, 27},
    {16126007, 28}, {16126008, 29}, {16126009, 30}};
}  // namespace StMyCuts

#endif