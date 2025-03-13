#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
//
// #include "HipoChain.h"
// #include "clas12reader.h"
//
// #include "src/classes/VetoHistograms/UpdateHistograms.cpp"
// #include "src/classes/clas12ana/clas12ana.cpp"
// #include "src/constants.h"
// #include "src/functions/Andrews_functions/Andrews_functions.cpp"
// #include "src/functions/GeneralFunctions.h"
// #include "src/functions/HipoChain_config.cpp"
// #include "src/functions/NeutronFunctions.h"
// #include "src/functions/neutron-veto/veto_functions.cpp"
#include "../../src/classes/HistPrinter/HistPrinter.cpp"
// #include "src/classes/HistPrinter/HistPrinter.cpp"

using namespace std;

void HistPrinterTester() {
    bool PrintOut = false;
    bool PrintOut1 = false;

    const char *filename = "/Users/alon/Downloads/Output_data_P9_run9_full/Andrew_plots_CD.root";
    // const char *filename = "/Users/alon/Downloads/Output_data_P9_run8_full_CutTester_test14/Andrew_plots_CD.root";
    TFile *file = new TFile(filename);

    string PDFFile = "/Users/alon/Downloads/TOut/TOut.pdf";
    system("rm -rf /Users/alon/Downloads/TOut");
    system("mkdir /Users/alon/Downloads/TOut");

    HistPrinter Printer;

    vector<TH1 *> HistoList;

    TString classname("TH1D");
    TString classnameTH2D("TH2D");

    TKey *Key;
    TIter Next((TList *)file->GetListOfKeys());

    while ((Key = (TKey *)Next())) {
        HistoList.push_back(((TH1 *)Key->ReadObj()));
        // if (Key->GetClassName() == classnameTH2D("TH1D")) {
        //     HistoList.push_back(((TH1D *) Key->ReadObj()));
        // } else if (Key->GetClassName() == classnameTH2D("TH2D")) {
        //     HistoList.push_back(((TH2D *) Key->ReadObj()));
        // }

        // if (Key->GetClassName() == classnameTH2D("TH1D")) {
        //     HistoList.push_back((TH1D *)Key);
        // } else if (Key->GetClassName() == classnameTH2D("TH2D")) {
        //     HistoList.push_back((TH2D *)Key);
        // }

        // string Histogram1DTempName = ((TH1D *)Key->ReadObj())->GetName();

        // if (PrintOut1) { cout << Histogram1DTempName << "\n\n"; }

        // if (findSubstring(Histogram1DTempName, Histogram1DNameSubstring) && (Key->GetClassName() != classnameTH2D("TH2D")) &&
        //     (Key->GetClassName() != classnameTFolder("TFolder")) && (Key->GetClassName() != classnameTHStack("THStack"))) {
        //     if (PrintOut) { cout << "\n\nKey name: " << ((TH1D *)Key->ReadObj())->GetName() << "; Type: " << Key->GetClassName() << "\n\n"; }

        //     string Histogram1DxLable = ((TH1D *)Key->ReadObj())->GetXaxis()->GetTitle();
        //     string Histogram1DTitle = ((TH1D *)Key->ReadObj())->GetTitle();

        //     if (PrintOut) {
        //         cout << "\nHistogram1DxLable = " << Histogram1DxLable << "\n";
        //         cout << "Histogram1DTitle = " << Histogram1DTitle << "\n";
        //         cout << "TLmom = " << TLmom << "\n";
        //     }

        //     if ((TLmom || !findSubstring(Histogram1DxLable, "Momentum"))) {
        //         HistogramFound = true;

        //         Histogram1D = ((TH1D *)Key->ReadObj());
        //         FoundHistName = Key->GetClassName();
        //         break;
        //     }
        // }
    }

    Printer.PlotHistograms(HistoList, PDFFile, false);
}