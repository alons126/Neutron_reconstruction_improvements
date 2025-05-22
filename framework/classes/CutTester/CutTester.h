//
// Created by Alon Sportes on 26/01/2025.
//

#ifndef CUTTESTER_H
#define CUTTESTER_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TTree.h"

// Include libraries:
#include "../../namespaces/general_utilities/utilities.h"

using namespace std;

class CutTester {
   private:
    TH1D *h_histo_BC_allN, *h_histo_BC_goodN, *h_histo_BC_badN;
    TH1D *h_histo_AC_allN, *h_histo_AC_goodN, *h_histo_AC_badN;

   public:
    // Constructor -----------------------------------------------------------------------------------------------------------------------------------------------------------

    CutTester() = default;

    // InitTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

    void InitTestHistograms(vector<TH1 *> &HistoList, const string &HistName, const string &HistTitle, const string &FinalState, const string &HistXLable, int NumberOfXBins, double LLim, double ULim);

    // FillTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

    void FillTestHistograms(bool isGN, bool isBN, double Variable, double weight, bool CutCondition);
};

#endif  // CUTTESTER_H
