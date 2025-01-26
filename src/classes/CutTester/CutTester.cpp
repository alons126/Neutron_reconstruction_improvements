//
// Created by Alon Sportes on 26/01/2025.
//

#include "CutTester.h"

// InitTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

void CutTester::InitTestHistograms(vector<TH1 *> &HistoList, const string &HistName, const string &HistTitle, const string &FinalState, const string &HistXLable, int NumberOfXBins,
                                   double LLim, double ULim) {
    string HistNameBC = HistName + "_BCTest_" + FinalState;
    string HistNameAC = HistName + "_ACTest_" + FinalState;
    string HistTitleBC = HistTitle + " Before Cut test;" + HistXLable;
    string HistTitleAC = HistTitle + " After Cut test;" + HistXLable;

    h_histo_BC = new TH1D(HistNameBC.c_str(), HistTitleBC.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_BC);
    h_histo_AC = new TH1D(HistNameAC.c_str(), HistTitleAC.c_str(), NumberOfXBins, LLim, ULim);
    HistoList.push_back(h_histo_AC);
}

// FillTestHistograms function -------------------------------------------------------------------------------------------------------------------------------------------

void CutTester::FillTestHistograms(double Variable, double weight, bool CutCondition) {
    h_histo_BC->Fill(Variable, weight);
    if (CutCondition) { h_histo_AC->Fill(Variable, weight); }
}
