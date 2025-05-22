//
// Created by Alon Sportes on 26/01/2025.
//

#include "CutTester.h"

// InitTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

void CutTester::InitTestHistograms(vector<TH1 *> &HistoList, const string &HistName, const string &HistTitle, const string &FinalState, const string &HistXLable, int NumberOfXBins, double LLim,
                                   double ULim) {
    string HistNameBC_allN = HistName + "_allN_BCTest_" + FinalState;
    string HistNameBC_goodN = HistName + "_goodN_BCTest_" + FinalState;
    string HistNameBC_badN = HistName + "_badN_BCTest_" + FinalState;

    string HistNameAC_allN = HistName + "_allN_ACTest_" + FinalState;
    string HistNameAC_goodN = HistName + "_goodN_ACTest_" + FinalState;
    string HistNameAC_badN = HistName + "_badN_ACTest_" + FinalState;

    string HistTitleBC_allN = HistTitle + " Before Cut Test (Good & Bad n);" + HistXLable + ";Counts";
    string HistTitleBC_goodN = HistTitle + " Before Cut Test (Good n Only);" + HistXLable + ";Counts";
    string HistTitleBC_badN = HistTitle + " Before Cut Test (Bad n Only);" + HistXLable + ";Counts";

    string HistTitleAC_allN = HistTitle + " After Cut Test (Good & Bad n);" + HistXLable + ";Counts";
    string HistTitleAC_goodN = HistTitle + " After Cut Test (Good n Only);" + HistXLable + ";Counts";
    string HistTitleAC_badN = HistTitle + " After Cut Test (Bad n Only);" + HistXLable + ";Counts";

    h_histo_BC_allN = new TH1D(HistNameBC_allN.c_str(), HistTitleBC_allN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_BC_allN);
    h_histo_AC_allN = new TH1D(HistNameAC_allN.c_str(), HistTitleAC_allN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_AC_allN);
    h_histo_BC_goodN = new TH1D(HistNameBC_goodN.c_str(), HistTitleBC_goodN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_BC_goodN);
    h_histo_AC_goodN = new TH1D(HistNameAC_goodN.c_str(), HistTitleAC_goodN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_AC_goodN);
    h_histo_BC_badN = new TH1D(HistNameBC_badN.c_str(), HistTitleBC_badN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_BC_badN);
    h_histo_AC_badN = new TH1D(HistNameAC_badN.c_str(), HistTitleAC_badN.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_AC_badN);
}

// FillTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

void CutTester::FillTestHistograms(bool isGN, bool isBN, double Variable, double weight, bool CutCondition) {
    h_histo_BC_allN->Fill(Variable, weight);

    if (isGN) {
        h_histo_BC_goodN->Fill(Variable, weight);
    } else if (isBN) {
        h_histo_BC_badN->Fill(Variable, weight);
    }

    if (CutCondition) {
        h_histo_AC_allN->Fill(Variable, weight);

        if (isGN) {
            h_histo_AC_goodN->Fill(Variable, weight);
        } else if (isBN) {
            h_histo_AC_badN->Fill(Variable, weight);
        }
    }
}
