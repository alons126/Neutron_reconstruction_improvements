//
// Created by Alon Sportes on 22/01/2025.
//

#ifndef HISTPRINTER_H
#define HISTPRINTER_H

#include <cstdlib>
#include <iostream>

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
#include "../../cuts/VetoCuts.h"

using namespace std;


class HistPrinter {
private:
    vector<const char *> summary_table_title = {"", "#(goodN)", "#(badN)", "#splitline{Signal}{Efficiency}", "#splitline{Signal}{Purity}"};

    vector<vector<const char *> > summary_table_bfSteps_epCDn;
    vector<vector<const char *> > summary_table_Step0_epCDn;
    vector<vector<const char *> > summary_table_Step1_epCDn;
    vector<vector<const char *> > summary_table_Step2_epCDn;

    vector<const char *> summary_table_bfSteps_epFDn;
    vector<const char *> summary_table_Step0_epFDn;
    vector<const char *> summary_table_Step1_epFDn;
    vector<const char *> summary_table_Step2_epFDn;

    vector<vector<const char *> > table_epCDn = {summary_table_title};
    bool First_table_epCDn_generation = true;

    vector<vector<const char *> > table_epFDn = {summary_table_title};

public:
    // Constructor
    // ======================================================================================================================================================================

    HistPrinter() = default;

    // PrintPage function
    // ======================================================================================================================================================================

    void PrintPage(const std::string &PageTitle, TCanvas *myText, char fileName[100], TLatex titles, TLatex text, const std::string &Constraint1,
                   const std::string &Constraint2);

    // GenerateSummaryTable function
    // ======================================================================================================================================================================

    void GenerateSummaryTable(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList,
                              TLatex titles, TLatex text, char fileName[100], string PDFFile, string Constraint1, string Constraint2,
                              bool LogScale2D);

    // SummaryTablePlotter function
    // ======================================================================================================================================================================

    void SummaryTablePlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList,
                             TLatex titles, TLatex text, char fileName[100], string PDFFile, string Constraint1, string Constraint2, bool LogScale2D);

    // GetHistogramEntries function
    // ======================================================================================================================================================================

    double GetHistogramEntries(const std::vector<TH1 *> HistoList, const std::string &histName);

    // extractStep function
    // ======================================================================================================================================================================

    std::string extractStep(const std::string &input);

    // SkippingCondition function
    // ======================================================================================================================================================================

    bool SkippingCondition(string HistoName, int canvas_ind);

    // replaceSubstring function
    // ======================================================================================================================================================================

    std::string replaceSubstring(const std::string &input, const std::string &toReplace, const std::string &replaceWith);

    // SectionPlotter function
    // ======================================================================================================================================================================

    void SectionPlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList,
                        string PDFFile, string Constraint1 = "", string Constraint2 = "", bool LogScale2D = false);

    // PlotHistograms function
    // ======================================================================================================================================================================

    void PlotHistograms(const vector<TH1 *> HistoList, const string &PDFFile, bool LogScale2D = false);
};


#endif //HISTPRINTER_H
