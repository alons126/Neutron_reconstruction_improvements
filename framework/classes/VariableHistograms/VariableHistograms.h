#ifndef VARIABLEHISTOGRAMS_H
#define VARIABLEHISTOGRAMS_H

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"

struct YVarConfig {
    std::string name;
    std::string title;
    std::string axisLabel;
    double min = 0.0;
    double max = 1.0;
    int bins = 50;
};

class VariableHistograms {
   public:
    enum NeutronQuality { All, Good, Bad };

    std::vector<TH1*> HistoList;
    std::set<std::string> skipList;
    std::string varName;

    TH1D* h_Var_allN = nullptr;
    TH1D* h_Var_goodN = nullptr;
    TH1D* h_Var_badN = nullptr;

    VariableHistograms(const std::string& varName_, const std::string& varTitle, const std::string& varAxisLabel, const std::string& finalStateTag, const std::vector<YVarConfig>& yVars)
        : varName(varName_) {
        make1D(varName_, varTitle, varAxisLabel, finalStateTag);

        for (const auto& y : yVars) { make2D(y, finalStateTag); }
    }

    ~VariableHistograms() {
        for (auto* h : HistoList) delete h;
    }

    void FillHistograms(double varVal, const std::map<std::string, double>& yvals, NeutronQuality nq) {
        // Choose suffix
        std::string suffix = (nq == All) ? "allN" : (nq == Good) ? "goodN" : "badN";

        // Fill 1D histogram
        std::string h1name = varName + "_" + suffix;
        for (auto* h : HistoList) {
            if (h->GetName() == h1name) {
                if (auto* h1 = dynamic_cast<TH1D*>(h)) h1->Fill(varVal);
                break;
            }
        }

        // Fill 2D histograms
        for (auto* h : HistoList) {
            auto* h2 = dynamic_cast<TH2D*>(h);
            if (!h2) continue;

            std::string name = h2->GetName();
            if (skipList.count(name)) continue;
            if (name.find("_" + suffix + "_") == std::string::npos) continue;

            size_t vsPos = name.find("_VS_");
            if (vsPos == std::string::npos) continue;
            std::string yvar = name.substr(2, vsPos - 2);  // skip "h_"

            auto it = yvals.find(yvar);
            if (it != yvals.end()) { h2->Fill(varVal, it->second); }
        }
    }

   private:
    void make1D(const std::string& var, const std::string& title, const std::string& label, const std::string& tag) {
        for (const auto& suffix : {"allN", "goodN", "badN"}) {
            std::string name = var + "_" + suffix + "_" + tag;
            if (gDirectory->Get(name.c_str())) {
                std::cerr << "Skipping 1D histogram: '" << name << "' already exists.\n";
                skipList.insert(name);
                continue;
            }

            std::string fullTitle = title + ";" + label + ";Counts";
            TH1D* h = new TH1D(name.c_str(), fullTitle.c_str(), 50, 0, 100);
            HistoList.push_back(h);

            if (suffix == std::string("allN"))
                h_Var_allN = h;
            else if (suffix == std::string("goodN"))
                h_Var_goodN = h;
            else
                h_Var_badN = h;
        }
    }

    void make2D(const YVarConfig& y, const std::string& tag) {
        for (const auto& suffix : {"allN", "goodN", "badN"}) {
            std::string name = "h_" + y.name + "_VS_" + varName + "_" + suffix + "_" + tag;

            if (gDirectory->Get(name.c_str())) {
                std::cerr << "Skipping 2D histogram: '" << name << "' already exists.\n";
                skipList.insert(name);
                continue;
            }

            std::string title = y.title + " vs " + varName + ";E^{CTOF}_{dep} [MeV];" + y.axisLabel;
            TH2D* h = new TH2D(name.c_str(), title.c_str(), 50, 0, 100, y.bins, y.min, y.max);
            HistoList.push_back(h);
        }
    }
};

#endif  // VARIABLEHISTOGRAMS_H