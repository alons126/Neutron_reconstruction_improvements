//
// Created by Alon Sportes on 13/03/2025.
//

#include "VariableHistograms.h"

using namespace std;
using namespace utilities;

// InitHistograms function -----------------------------------------------------------------------------------------------------------------------------------------------

void VariableHistograms::InitHistograms(vector<vector<TH1 *>> &HistoList, const string VarName, const string Step, const string FinalState, const string TitleAdditions, const bool separateGoodBadN) {
    OnlySeparatedN = separateGoodBadN;

    vector<TH1 *> TempHistoList;

    string NameEnding = (Step != "") ? "_" + Step + "_" + FinalState : "_" + FinalState;
    string TitleEnding = (TitleAdditions != "") ? " (" + TitleAdditions + ")" : "";

    if (OnlySeparatedN) {
        // 1D histogram
        h_Var_1D = new TH1D((VariableNames[VarName]["VarName"] + NameEnding).c_str(), (VariableNames[VarName]["VarLabel"] + " distribution" + TitleAdditions).c_str(), 50, 0, 1);
        TempHistoList.push_back(h_Var_1D);

        // 2D histograms
        h_P_n_VS_Var_2D = new TH2D(("P_n_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                   (VariableNames["P_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_n_VS_Var_2D);

        h_theta_n_VS_Var_2D = new TH2D(("theta_n_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                       (VariableNames["theta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_VS_Var_2D);

        h_phi_n_VS_Var_2D = new TH2D(("phi_n_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                     (VariableNames["phi_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_n_VS_Var_2D);

        h_P_miss_VS_Var_2D = new TH2D(("P_miss_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                      (VariableNames["P_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_miss_VS_Var_2D);

        h_theta_miss_VS_Var_2D = new TH2D(("theta_miss_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                          (VariableNames["theta_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_miss_VS_Var_2D);

        h_phi_miss_VS_Var_2D = new TH2D(("phi_miss_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                        (VariableNames["phi_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_miss_VS_Var_2D);

        h_dpp_VS_Var_2D = new TH2D(("dpp_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                   (VariableNames["dpp"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_dpp_VS_Var_2D);

        h_theta_n_miss_VS_Var_2D = new TH2D(("theta_n_miss_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                            (VariableNames["theta_n_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_miss_VS_Var_2D);

        h_beta_n_VS_Var_2D = new TH2D(("beta_n_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                      (VariableNames["beta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_beta_n_VS_Var_2D);

        h_E_p_VS_Var_2D = new TH2D(("E_p_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                   (VariableNames["E_p"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_p_VS_Var_2D);

        h_E_miss_VS_Var_2D = new TH2D(("E_miss_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                      (VariableNames["E_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_miss_VS_Var_2D);

        h_path_VS_Var_2D = new TH2D(("path_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                    (VariableNames["path"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_path_VS_Var_2D);

        h_ToF_VS_Var_2D = new TH2D(("ToF_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                   (VariableNames["ToF"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_ToF_VS_Var_2D);

        h_nSector_VS_Var_2D = new TH2D(("nSector_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                       (VariableNames["nSector"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_nSector_VS_Var_2D);

        h_Edep_CND_VS_Var_2D = new TH2D(("Edep_CND_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                        (VariableNames["Edep_CND"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND_VS_Var_2D);

        h_Edep_CND1_VS_Var_2D = new TH2D(("Edep_CND1_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Edep_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND1_VS_Var_2D);

        h_Edep_CND2_VS_Var_2D = new TH2D(("Edep_CND2_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Edep_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND2_VS_Var_2D);

        h_Edep_CND3_VS_Var_2D = new TH2D(("Edep_CND3_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Edep_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND3_VS_Var_2D);

        h_Size_CND1_VS_Var_2D = new TH2D(("Size_CND1_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Size_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND1_VS_Var_2D);

        h_Size_CND2_VS_Var_2D = new TH2D(("Size_CND2_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Size_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND2_VS_Var_2D);

        h_Size_CND3_VS_Var_2D = new TH2D(("Size_CND3_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                         (VariableNames["Size_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND3_VS_Var_2D);

        h_LayerMulti_CND1_VS_Var_2D = new TH2D(("LayerMulti_CND1_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                               (VariableNames["LayerMulti_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND1_VS_Var_2D);

        h_LayerMulti_CND2_VS_Var_2D = new TH2D(("LayerMulti_CND2_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                               (VariableNames["LayerMulti_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND2_VS_Var_2D);

        h_LayerMulti_CND3_VS_Var_2D = new TH2D(("LayerMulti_CND3_VS_" + VariableNames[VarName]["VarName"] + NameEnding).c_str(),
                                               (VariableNames["LayerMulti_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND3_VS_Var_2D);
    } else {
        // 1D histogram
        h_Var_goodN_1D = new TH1D((VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(), (VariableNames[VarName]["VarLabel"] + " distribution" + TitleAdditions).c_str(), 50, 0, 1);
        TempHistoList.push_back(h_Var_goodN_1D);

        h_Var_badN_1D = new TH1D((VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(), (VariableNames[VarName]["VarLabel"] + " distribution" + TitleAdditions).c_str(), 50, 0, 1);
        TempHistoList.push_back(h_Var_badN_1D);

        // 2D histograms
        h_P_n_VS_Var_goodN_2D = new TH2D(("P_n_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                         (VariableNames["P_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_n_VS_Var_goodN_2D);

        h_P_n_VS_Var_badN_2D = new TH2D(("P_n_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                        (VariableNames["P_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_n_VS_Var_badN_2D);

        h_theta_n_VS_Var_goodN_2D = new TH2D(("theta_n_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                             (VariableNames["theta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_VS_Var_goodN_2D);

        h_theta_n_VS_Var_badN_2D = new TH2D(("theta_n_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                            (VariableNames["theta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_VS_Var_badN_2D);

        h_phi_n_VS_Var_goodN_2D = new TH2D(("phi_n_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                           (VariableNames["phi_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_n_VS_Var_goodN_2D);

        h_phi_n_VS_Var_badN_2D = new TH2D(("phi_n_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                          (VariableNames["phi_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_n_VS_Var_badN_2D);

        h_P_miss_VS_Var_goodN_2D = new TH2D(("P_miss_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                            (VariableNames["P_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_miss_VS_Var_goodN_2D);

        h_P_miss_VS_Var_badN_2D = new TH2D(("P_miss_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                           (VariableNames["P_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_P_miss_VS_Var_badN_2D);

        h_theta_miss_VS_Var_goodN_2D = new TH2D(("theta_miss_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                                (VariableNames["theta_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_miss_VS_Var_goodN_2D);

        h_theta_miss_VS_Var_badN_2D = new TH2D(("theta_miss_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                               (VariableNames["theta_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_miss_VS_Var_badN_2D);

        h_phi_miss_VS_Var_goodN_2D = new TH2D(("phi_miss_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                              (VariableNames["phi_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_miss_VS_Var_goodN_2D);

        h_phi_miss_VS_Var_badN_2D = new TH2D(("phi_miss_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                             (VariableNames["phi_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_phi_miss_VS_Var_badN_2D);

        h_dpp_VS_Var_goodN_2D = new TH2D(("dpp_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                         (VariableNames["dpp"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_dpp_VS_Var_goodN_2D);

        h_dpp_VS_Var_badN_2D = new TH2D(("dpp_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                        (VariableNames["dpp"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_dpp_VS_Var_badN_2D);

        h_theta_n_miss_VS_Var_goodN_2D = new TH2D(("theta_n_miss_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                                  (VariableNames["theta_n_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_miss_VS_Var_goodN_2D);

        h_theta_n_miss_VS_Var_badN_2D = new TH2D(("theta_n_miss_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                                 (VariableNames["theta_n_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_theta_n_miss_VS_Var_badN_2D);

        h_beta_n_VS_Var_goodN_2D = new TH2D(("beta_n_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                            (VariableNames["beta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_beta_n_VS_Var_goodN_2D);

        h_beta_n_VS_Var_badN_2D = new TH2D(("beta_n_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                           (VariableNames["beta_n"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_beta_n_VS_Var_badN_2D);

        h_E_p_VS_Var_goodN_2D = new TH2D(("E_p_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                         (VariableNames["E_p"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_p_VS_Var_goodN_2D);

        h_E_p_VS_Var_badN_2D = new TH2D(("E_p_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                        (VariableNames["E_p"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_p_VS_Var_badN_2D);

        h_E_miss_VS_Var_goodN_2D = new TH2D(("E_miss_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                            (VariableNames["E_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_miss_VS_Var_goodN_2D);

        h_E_miss_VS_Var_badN_2D = new TH2D(("E_miss_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                           (VariableNames["E_miss"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_E_miss_VS_Var_badN_2D);

        h_path_VS_Var_goodN_2D = new TH2D(("path_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                          (VariableNames["path"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_path_VS_Var_goodN_2D);

        h_path_VS_Var_badN_2D = new TH2D(("path_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                         (VariableNames["path"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_path_VS_Var_badN_2D);

        h_ToF_VS_Var_goodN_2D = new TH2D(("ToF_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                         (VariableNames["ToF"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_ToF_VS_Var_goodN_2D);

        h_ToF_VS_Var_badN_2D = new TH2D(("ToF_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                        (VariableNames["ToF"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_ToF_VS_Var_badN_2D);

        h_nSector_VS_Var_goodN_2D = new TH2D(("nSector_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                             (VariableNames["nSector"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_nSector_VS_Var_goodN_2D);

        h_nSector_VS_Var_badN_2D = new TH2D(("nSector_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                            (VariableNames["nSector"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_nSector_VS_Var_badN_2D);

        h_Edep_CND_VS_Var_goodN_2D = new TH2D(("Edep_CND_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                              (VariableNames["Edep_CND"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND_VS_Var_goodN_2D);

        h_Edep_CND_VS_Var_badN_2D = new TH2D(("Edep_CND_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                             (VariableNames["Edep_CND"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND_VS_Var_badN_2D);

        h_Edep_CND1_VS_Var_goodN_2D = new TH2D(("Edep_CND1_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Edep_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND1_VS_Var_goodN_2D);

        h_Edep_CND1_VS_Var_badN_2D = new TH2D(("Edep_CND1_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Edep_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND1_VS_Var_badN_2D);

        h_Edep_CND2_VS_Var_goodN_2D = new TH2D(("Edep_CND2_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Edep_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND2_VS_Var_goodN_2D);

        h_Edep_CND2_VS_Var_badN_2D = new TH2D(("Edep_CND2_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Edep_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND2_VS_Var_badN_2D);

        h_Edep_CND3_VS_Var_goodN_2D = new TH2D(("Edep_CND3_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Edep_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND3_VS_Var_goodN_2D);

        h_Edep_CND3_VS_Var_badN_2D = new TH2D(("Edep_CND3_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Edep_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Edep_CND3_VS_Var_badN_2D);

        h_Size_CND1_VS_Var_goodN_2D = new TH2D(("Size_CND1_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Size_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND1_VS_Var_goodN_2D);

        h_Size_CND1_VS_Var_badN_2D = new TH2D(("Size_CND1_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Size_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND1_VS_Var_badN_2D);

        h_Size_CND2_VS_Var_goodN_2D = new TH2D(("Size_CND2_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Size_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND2_VS_Var_goodN_2D);

        h_Size_CND2_VS_Var_badN_2D = new TH2D(("Size_CND2_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Size_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND2_VS_Var_badN_2D);

        h_Size_CND3_VS_Var_goodN_2D = new TH2D(("Size_CND3_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                               (VariableNames["Size_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND3_VS_Var_goodN_2D);

        h_Size_CND3_VS_Var_badN_2D = new TH2D(("Size_CND3_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                              (VariableNames["Size_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_Size_CND3_VS_Var_badN_2D);

        h_LayerMulti_CND1_VS_Var_goodN_2D = new TH2D(("LayerMulti_CND1_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                                     (VariableNames["LayerMulti_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND1_VS_Var_goodN_2D);

        h_LayerMulti_CND1_VS_Var_badN_2D = new TH2D(("LayerMulti_CND1_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                                    (VariableNames["LayerMulti_CND1"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND1_VS_Var_badN_2D);

        h_LayerMulti_CND2_VS_Var_goodN_2D = new TH2D(("LayerMulti_CND2_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                                     (VariableNames["LayerMulti_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND2_VS_Var_goodN_2D);

        h_LayerMulti_CND2_VS_Var_badN_2D = new TH2D(("LayerMulti_CND2_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                                    (VariableNames["LayerMulti_CND2"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND2_VS_Var_badN_2D);

        h_LayerMulti_CND3_VS_Var_goodN_2D = new TH2D(("LayerMulti_CND3_VS_" + VariableNames[VarName]["VarName"] + "_goodN" + NameEnding).c_str(),
                                                     (VariableNames["LayerMulti_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND3_VS_Var_goodN_2D);

        h_LayerMulti_CND3_VS_Var_badN_2D = new TH2D(("LayerMulti_CND3_VS_" + VariableNames[VarName]["VarName"] + "_badN" + NameEnding).c_str(),
                                                    (VariableNames["LayerMulti_CND3"]["VarLabel"] + " vs. " + VariableNames[VarName]["VarLabel"] + TitleAdditions).c_str(), 100, 0, 1, 100, 0, 1);
        TempHistoList.push_back(h_LayerMulti_CND3_VS_Var_badN_2D);
    }

    HistoList.push_back(TempHistoList);
}
