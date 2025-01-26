//
// Created by Alon Sportes on 26/01/2025.
//

#ifndef CUTTESTER_H
#define CUTTESTER_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../../functions/GeneralFunctions.h"
#include "../../constants.h"

using namespace std;

class CutTester {
private:
    TH1D *h_histo_BC;
    TH1D *h_histo_AC;

public:
    // Constructor ----------------------------------------------------------------------------------------------------------------------------------------------------------

    CutTester() = default;

    // InitTestHistograms function ------------------------------------------------------------------------------------------------------------------------------------------

    void InitTestHistograms(vector<TH1 *> HistoList, const string &HistName, const string &HistTitle, const string &FinalState, const string &HistXLable, int NumberOfXBins, double LLim,
                                   double ULim);

    // FillTestHistograms function ------------------------------------------------------------------------------------------------------------------------------------------

    void FillTestHistograms(double Variable, double weight, bool CutCondition);
};


#endif //CUTTESTER_H
