#ifndef StMyCuts_h
#define StMyCuts_h

#include "Rtypes.h"
#include <map>
#include <vector>

namespace StMyCuts{
	const std::vector<UInt_t> Trigger   = {720007};
	const std::vector<Double_t> VzCut   = {198., 202.};
	const Double_t VrCut                = 2.;
	const std::vector<Double_t> VrCen   = {0., -2.};
	const Double_t DcaCut               = 3.;
	const Int_t    NHitsFitCut          = 10;
	const Double_t NHitsRatioCut        = 0.52;
	const Int_t    NHitsDedxCut         = 5;
	const Double_t NSigmaCut            = 2.;
	const std::vector<Double_t> PtCutQ  = {0.2, 2.};
	const std::vector<Double_t> EtaCutQ = {-2., -1.1, -1., 0.};

	const std::map<Int_t, Int_t> RunIdIndex = {{20355020, 0}, {20355021, 1}, {21044023, 2}, {21044024, 3}, {21044025, 4}, {21044027, 5}, {21044028, 6}, {21044029, 7}, {21044030, 8}, {21044031, 9}, {21044032, 10}, {21044033, 11}, {21044034, 12}, {21044036, 13}, {21044037, 14}, {21044038, 15}, {21044039, 16}, {21044040, 17}, {21044041, 18}, {21044042, 19}, {21044043, 20}, {21045001, 21}, {21045002, 22}, {21045003, 23}, {21045005, 24}, {21045006, 25}, {21045007, 26}, {21045008, 27}, {21045009, 28}, {21045010, 29}, {21045011, 30}};
}

#endif