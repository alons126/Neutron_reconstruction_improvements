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

void HipoChain_config(HipoChain &chain, const string &sn, const string &AnalyseFilePath, const string &AnalyseFileSample, const string &AnalyseFile)
{
    /* Data in cache/clas12/rg-m/production/pass1/2gev/D/dst/recon */
    vector<string> Runs = {
        "015567", "015573", "015578", "015583", "015590", "015595", "015602", "015608", "015613", "015618", "015624",
        "015568", "015574", "015579", "015586", "015591", "015598", "015603", "015609", "015614", "015619", "015625",
        "015569", "015575", "015580", "015587", "015592", "015599", "015604", "015610", "015615", "015620", "015626",
        "015570", "015576", "015581", "015588", "015593", "015600", "015606", "015611", "015616", "015622", "015627",
        "015572", "015577", "015582", "015589", "015594", "015601", "015607", "015612", "015617", "015623"};

    for (int i = 0; i < Runs.size(); i++)
    {
        string TempAnalyseFile = "/" + AnalyseFilePath + "/" + Runs.at(i) + "/*.hipo";
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
}