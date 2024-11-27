#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "HipoChain.h"
// #include "eventcut.h"
// #include "functions.h"

using namespace std;
using namespace clas12;

void HipoChain_config(HipoChain &chain, const string &AnalyseFilePath)
{
    if (AnalyseFilePath == "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*")
    {
        const bool PrintOut = true;

        string D2_6GeV_Data_Path = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/";

        /* Data in cache/clas12/rg-m/production/pass1/6gev/D/dst/recon */
        vector<string> Runs = {
            "015045", "015052", "015058", "015066", "015077", "015094", "015100", "015106", "015442", "015449", "015456",
            "015046", "015053", "015059", "015067", "015078", "015095", "015101", "015435", "015443", "015450",
            "015047", "015054", "015060", "015072", "015079", "015096", "015102", "015436", "015444", "015451",
            "015049", "015055", "015061", "015073", "015081", "015097", "015103", "015437", "015445", "015452",
            "015050", "015056", "015062", "015074", "015082", "015098", "015104", "015439", "015447", "015454",
            "015051", "015057", "015065", "015075", "015093", "015099", "015105", "015441", "015448", "015455"};

        for (int i = 0; i < Runs.size(); i++)
        {
            string TempAnalyseFile = D2_6GeV_Data_Path + Runs.at(i) + "/*.hipo";
            chain.Add(TempAnalyseFile.c_str());

            if (PrintOut)
            {
                cout << TempAnalyseFile << " directory added to HipoChain!\n";
            }
        }

        if (PrintOut)
        {
            cout << "\n";
        }
    } else {
        chain.Add(AnalyseFilePath);
    }
}